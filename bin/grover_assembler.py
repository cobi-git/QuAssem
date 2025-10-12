#!/usr/bin/env python3

# --------------------------------------
#  Author: Inuk Jung, inukjung@knu.ac.kr (Kyungpook National University)
#  Co-author: Yongtae Kim, yongtae@knu.ac.kr (Kyungpook National University)
#  Date: 2025-10-11
#
#  Copyright (c) 2025 Inuk Jung & Yongtae Kim
#  All rights reserved.
# --------------------------------------

# The Grover-DP assembler
#
# How to run:
#   $ python grover_assembler.py
#
# Note: 
# - Run the program again if the assembled result is incorrect
# - This can happen because the k-mer size is small and the probabilistics are 
#   not strong enough even after the diffusing process
# - Setting the GROVER_SHOTS or the max_resamples higher may help but not by much
# - Will improve this Grover-DP in the upcoming version by increasing the k-mer size 
#   and making the entire algorithm more robust

from pathlib import Path
import math

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator

# Paths
SCRIPT_PATH = Path(__file__).resolve()
BIN_DIR = SCRIPT_PATH.parent
PROJECT_DIR = BIN_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RES_DIR = PROJECT_DIR / "res"

# Input and output paths
IN_FASTA   = DATA_DIR / "sarscov_N_4mer_shuffled.fasta"
OUT_FASTA  = RES_DIR  / "sarscov_N_4mer_shuffled_assembled.fasta"

# parmeters used for assembly and the quantum circuit
# Note: setting k-mer size to 4 for small qubits and easier drawing of the circuit diagram
KMER = 4
KMER_BITS = 2 * KMER
GROVER_SHOTS = 16   # no. of repeating the circuit
LINE_WRAP = 70

# Grover object
class GroverProgress:
    def __init__(self):
        self.calls = 0
    
# Reading sequences from FASTA file
def read_fasta(path):
    p = Path(path)
    recs = []
    cur_id = None
    cur_seq = []
    with p.open("r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    recs.append((cur_id, "".join(cur_seq).upper()))
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        
        if cur_id is not None:
            recs.append((cur_id, "".join(cur_seq).upper()))

    print(f"[loaded] {str(p)} file")
    return recs

# writing final assembly result in FASTA format
def write_fasta(records, out_path, line_wrap):
    out_path = Path(out_path)
    with out_path.open("w") as fw:
        for header, seq in records:
            fw.write(f">{header}\n")
            if line_wrap and line_wrap > 0:
                for i in range(0, len(seq), line_wrap):
                    fw.write(seq[i:i+line_wrap] + "\n")
            else:
                fw.write(seq + "\n")

    return out_path

# printing assembled sequences in FASTA format
def print_fasta(records, line_wrap):
    for header, seq in records:
        print(f">{header}")
        if line_wrap and line_wrap > 0:
            for i in range(0, len(seq), line_wrap):
                print(seq[i:i+line_wrap])
        else:
            print(seq)


# 2-bit encoding scheme
# Note: A-00, C-01, G-10, T-11
NT2BITS = {"A": (0,0), "C": (0,1), "G": (1,0), "T": (1,1)}

# encoding kmer to 2-bit encodings
def encode_kmer_bits(s):
    out = []
    for ch in s:
        b0, b1 = NT2BITS[ch]
        out.extend([b0, b1])  # per-base 2 bits
    return out

# encoding prefix k-mer sequence
def prefix_bits(seq):
    return encode_kmer_bits(seq[:KMER])

# encoding suffix k-mer sequence
def suffix_bits(seq):
    return encode_kmer_bits(seq[-KMER:])

# converting an integer x to n little-endian bits
# Note: used by the QROM for addressing the sequence indexes
def int_to_bits_le(x, n):
    return [(x >> i) & 1 for i in range(n)]

# apply the MCX to flip target iff all qubits in controls are |1⟩
# Note: used by QROM loader, oracle and diffuser
def _mcx_noancilla(circ, controls, target):
    circ.mcx(controls, target)

# QROM loader (index-controlled) functions
def apply_index_controlled_x(circ, idx_q, tgt, index_val):
    n = len(idx_q)
    bits = int_to_bits_le(index_val, n)

    for qb, bit in zip(idx_q, bits):
        if bit == 0:
            circ.x(qb)
    if n == 0:
        circ.x(tgt)
    else:
        _mcx_noancilla(circ, list(idx_q), tgt)
    for qb, bit in zip(idx_q, bits):
        if bit == 0:
            circ.x(qb)

def make_qrom_loader(active_kmers, n_idx):
    k = KMER_BITS
    idx_q  = QuantumRegister(n_idx, "idx")
    kmer_q = QuantumRegister(k, "kmer")
    sub = QuantumCircuit(idx_q, kmer_q, name="QROM_LOAD")

    M = len(active_kmers)
    for i in range(M):
        bits = active_kmers[i]
        for b in range(k):
            if bits[b] == 1:
                apply_index_controlled_x(sub, idx_q, kmer_q[b], i)
    return sub

# Oracle component
def apply_equality_phase_oracle(circ, kmer_q, oracle_q, target_bits):
    for q, t in zip(kmer_q, target_bits):
        if t == 0:
            circ.x(q)
    controls = list(kmer_q)
    if len(controls) == 0:
        circ.x(oracle_q)
    else:
        _mcx_noancilla(circ, controls, oracle_q)
    for q, t in zip(kmer_q, target_bits):
        if t == 0:
            circ.x(q)

# Diffuser component
def apply_diffuser(circ, idx_q):
    n = len(idx_q)
    circ.h(idx_q)
    circ.x(idx_q)
    if n == 1:
        circ.z(idx_q[0])
    else:
        circ.h(idx_q[-1])
        _mcx_noancilla(circ, list(idx_q[:-1]), idx_q[-1])
        circ.h(idx_q[-1])
    circ.x(idx_q)
    circ.h(idx_q)

# building the grover circuit
def build_grover_circuit(active_kmers, target_bits, num_iters):
    M = len(active_kmers)
    n_idx = max(1, math.ceil(math.log2(max(1, M))))
    idx_q   = QuantumRegister(n_idx, "idx")
    kmer_q  = QuantumRegister(KMER_BITS, "kmer")
    oracle  = QuantumRegister(1, "ora")
    creg    = ClassicalRegister(n_idx, "cidx")
    qc = QuantumCircuit(idx_q, kmer_q, oracle, creg, name=f"grover_M{M}_R{num_iters}")

    # uniform over indices + prepare |-> on oracle
    qc.h(idx_q)
    qc.x(oracle); qc.h(oracle)

    # build QROM loader as a subcircuit (NOT converted to a Gate)
    loader_circ = make_qrom_loader(active_kmers, n_idx)
    loader_inv  = loader_circ.inverse()
    loader_qubits = list(idx_q) + list(kmer_q)

    R = max(1, num_iters)
    for _ in range(R):
        # inline the loader (no custom instruction emitted)
        qc.compose(loader_circ, loader_qubits, inplace=True)

        # oracle: phase flip iff loaded kmer == target
        apply_equality_phase_oracle(qc, kmer_q, oracle[0], target_bits)

        # uncompute loader (inline inverse)
        qc.compose(loader_inv, loader_qubits, inplace=True)

        # diffuser on index
        apply_diffuser(qc, idx_q)

    # unprepare oracle and measure
    qc.h(oracle); qc.x(oracle)
    qc.measure(idx_q, creg)
    return qc


# The main grover algorithm (the extension process based on the k-mer matching pairs)
# Note: returns a measured integer in [0, M), or None on failure.
# Note: equivalent to the @find_best_extension_for function in the classical assembler
def grover_pick_index(sim, active_kmers, target_bits, *, label, shots=GROVER_SHOTS, max_resamples=20):
    M = len(active_kmers)
    if M == 0:
        print(f"[grover:{label}] no active candidates")
        return None

    n_idx = max(1, math.ceil(math.log2(max(1, M))))
    R = max(1, int(math.floor(math.pi/4 * math.sqrt(M))))

    for attempt in range(1, max_resamples + 1):
        print(f"[grover:{label}] attempt {attempt}/{max_resamples} M={M} n_idx={n_idx} "
              f"kmer_bits={KMER_BITS} R={R} shots={shots}", flush=True)

        qc = build_grover_circuit(active_kmers, target_bits, num_iters=R)
        res = sim.run(qc, shots=shots).result()

        try:
            ops = qc.count_ops()
            ops_info = ", ".join(f"{k}:{v}" for k, v in list(ops.items())[:6])
        except Exception:
            ops_info = "na"
        print(f"[grover:{label}] {ops_info}", flush=True)

        counts = res.get_counts()
        if not counts:
            print(f"[grover:{label}] no counts returned")
            continue

        winner, wcount = max(counts.items(), key=lambda kv: kv[1])
        idx_val = int(winner, 2)  # IMPORTANT: no reverse
        print(f"[grover:{label}] winner={winner} ({wcount}/{shots}) → idx_rel={idx_val}", flush=True)

        if 0 <= idx_val < M:
            return idx_val

        print(f"[grover:{label}] measured out-of-range; retrying...", flush=True)

    print(f"[grover:{label}] exhausted {max_resamples} attempts without a valid index", flush=True)
    return None

# object for storing alignment results
class Contig:
    __slots__ = ("seq", "ids")
    def __init__(self, seq, ids):
        self.seq = seq
        self.ids = list(ids)

# Exact K-mer comparison by the suffix (hashing, so O(1) time)
def exact_kmer_right(a, b, k):
    return len(a) >= k and len(b) >= k and a[-k:] == b[:k]

# Exact K-mer comparison by the prefix (hashing, so O(1) time)
def exact_kmer_left(a, b, k):
    return len(a) >= k and len(b) >= k and b[-k:] == a[:k]

# the main assembly routine
# merging sequences to build the assembled contig
def assemble_with_grover(reads):
    ids  = [r[0] for r in reads]
    seqs = [r[1] for r in reads]
    N = len(seqs)

    pref_bits = [prefix_bits(s) for s in seqs]
    suff_bits = [suffix_bits(s) for s in seqs]

    remaining = set(range(N))
    contigs = []

    sim = AerSimulator(method="statevector")

    round_no = 0
    while remaining:
        round_no += 1
        i0 = min(remaining)
        contig = Contig(seqs[i0], [ids[i0]])
        remaining.remove(i0)
        print(f"[seed] {ids[i0]} | remaining={len(remaining)}", flush=True)

        extended = True
        while extended and remaining:
            extended = False
            active = sorted(list(remaining))
            M = len(active)

            print(f"[status] round={round_no} contig_len={len(contig.seq)} "
                  f"reads_used={N-len(remaining)}/{N} active={M}", flush=True)

            # RIGHT extension (suffix(contig) vs prefix(active))
            target = suffix_bits(contig.seq)
            active_kmers = [pref_bits[j] for j in active]
            print(f"[try] RIGHT: target=suffix({contig.ids[-1]}) vs prefix(active)", flush=True)
            idx_rel = grover_pick_index(sim, active_kmers, target_bits=target, label="RIGHT")

            if idx_rel is not None:
                j = active[idx_rel]
                if exact_kmer_right(contig.seq, seqs[j], KMER):
                    contig.seq = contig.seq + seqs[j][KMER:]
                    contig.ids.append(ids[j])
                    remaining.remove(j)
                    print(f"[merge] RIGHT {contig.ids[-2]} -> {ids[j]} | new_len={len(contig.seq)} "
                          f"| remaining={len(remaining)}", flush=True)
                    extended = True
                else:
                    print("[right] Grover pick failed verification", flush=True)
            else:
                print("[right] no valid RIGHT index from Grover", flush=True)

            if extended or not remaining:
                continue  # if we extended RIGHT (or nothing left), skip LEFT this round

            # LEFT extension (suffix(active) vs prefix(contig))
            active = sorted(list(remaining))
            M = len(active)
            target = prefix_bits(contig.seq)
            active_kmers = [suff_bits[jj] for jj in active]
            print(f"[try] LEFT : target=prefix({contig.ids[0]}) vs suffix(active)", flush=True)
            idx_rel = grover_pick_index(sim, active_kmers, target_bits=target, label="LEFT")

            if idx_rel is not None:
                j = active[idx_rel]
                if exact_kmer_left(contig.seq, seqs[j], KMER):
                    contig.seq = seqs[j] + contig.seq[KMER:]
                    contig.ids = [ids[j]] + contig.ids
                    remaining.remove(j)
                    print(f"[merge] LEFT {ids[j]} -> {contig.ids[1]} | new_len={len(contig.seq)} "
                          f"| remaining={len(remaining)}", flush=True)
                    extended = True
                else:
                    print("[left] Grover pick failed verification", flush=True)
            else:
                print("[left] no valid LEFT index from Grover", flush=True)

        print(f"[finalized contig] len={len(contig.seq)} reads={len(contig.ids)} "
              f"order={'->'.join(contig.ids)}", flush=True)
        contigs.append(contig)

    print(f"[done] total contigs={len(contigs)}", flush=True)
    return contigs

# setting the contid id and information of its assembly (includes order of assembled seqs)
def contigs_to_fasta_records(contigs):
    recs = []
    for i, c in enumerate(contigs, start=1):
        header = f"contig_{i} nreads={len(c.ids)} len={len(c.seq)} order={'->'.join(c.ids)}"
        recs.append((header, c.seq))
    return recs

# main execution
def main():
    print("=== Grover-DP Assembler ===")
    print(f"[project] {PROJECT_DIR}")
    print(f"[data]    {DATA_DIR}")
    print(f"[input]   {IN_FASTA}")
    print(f"[res]     {RES_DIR}")

    # read input FASTA file
    reads = read_fasta(IN_FASTA)
    print(f"[loaded] {len(reads)} reads")

    # perform the assembly
    print("\n--- Assembly START ---")
    contigs = assemble_with_grover(reads)
    print("\n--- Assembly END ---")

    # store assembed result in FASTA format
    fasta_records = contigs_to_fasta_records(contigs)

    # printing assembly results to stdout
    print("\n--- Assembled Contig ---")
    print_fasta(fasta_records, line_wrap=LINE_WRAP)

    # write assembled result to file
    out_fasta = RES_DIR / "sarscov_S_4mer_seqs_shuffled_grover_assembled.fasta"
    write_fasta(fasta_records, out_fasta, line_wrap=LINE_WRAP)

    print("\n--- Assembly summary and order ---")
    print(f"[contigs] assembled contig written to {out_fasta}")
    for i, c in enumerate(contigs, start=1):
        print(f"  - contig_{i}: len={len(c.seq)} reads={len(c.ids)}")
        print(f"  - order={'->'.join(c.ids)}")

if __name__ == "__main__":
    main()