#!/usr/bin/env python3

# --------------------------------------
#  Author: Inuk Jung, inukjung@knu.ac.kr (Kyungpook National University)
#  Co-author: Yongtae Kim, yongtae@knu.ac.kr (Kyungpook National University)
#  Date: 2025-10-11
#
#  Copyright (c) 2025 Inuk Jung & Yongtae Kim
#  All rights reserved.
# --------------------------------------

# A simple classical DP based assembler
#
# How to run:
#   $ python classical_assembler.py
#
# Note:
# - This assembler is a constrained version of the Smith-Waterman algorithm 
#   specialized for the extension on the suffix/prefix matching (overlapping) sequence pairs

from pathlib import Path
from dataclasses import dataclass

# Paths
SCRIPT_PATH = Path(__file__).resolve()
BIN_DIR     = SCRIPT_PATH.parent
PROJECT_DIR = BIN_DIR.parent
DATA_DIR    = PROJECT_DIR / "data"
RES_DIR     = PROJECT_DIR / "res"

# Input and output paths
IN_FASTA    = DATA_DIR / "sarscov_N_10mer_shuffled.fasta"
OUT_FASTA   = RES_DIR / "sarscov_N_10mer_shuffled_assembled.fasta"

# parmeters used for assembly
KMER = 10
MATCH_SCORE = 1
MISMATCH_PENALTY = -1
GAP_PENALTY = -1
LINE_WRAP = 70
VERBOSE = True

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

# object for storing alignment results
@dataclass(frozen=True)
class AlignmentRecord:
    score: int
    i_start: int
    i_end:   int
    j_start: int
    j_end:   int
    aln1:    str
    aln2:    str

# The DP alignment function
# A: the full length query sequence. B: the full length target sequence (was checked by the CMX gate in the Grover-DP assembler)
# Note: Sequence A and B here are not K-mer sequences
def overlap_alignment_suffix_prefix(A, B):
    n, m = len(A), len(B)
    NEG_INF = -10**9
    H = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    P = [[0] * (m + 1) for _ in range(n + 1)]  # 0 none, 1 diag (match), 2 up (mismatch), 3 left (mismatch)

    # initializing the scoring matrix
    H[0][0] = 0 # the left top cell (staring point for DP scoring)

    # initializing the rows
    for j in range(1, m + 1):
        H[0][j] = NEG_INF   # instead of 0, we set the top row with NEG_INF or -∞ so that it doesnt skip characters at the sart of the B sequences

    # initializing the columns
    for i in range(1, n + 1):
        H[i][0] = 0

    # performing the DP scoring (classical scoring SW scoring scheme)
    for i in range(1, n + 1):
        ai = A[i - 1]
        for j in range(1, m + 1):
            bj = B[j - 1]
            diag = H[i - 1][j - 1] + (MATCH_SCORE if ai == bj else MISMATCH_PENALTY) # match
            up   = H[i - 1][j] + GAP_PENALTY # mismatch
            left = H[i][j - 1] + GAP_PENALTY # mismatch
            if diag >= up and diag >= left:
                H[i][j] = diag; P[i][j] = 1
            elif up >= left:
                H[i][j] = up;   P[i][j] = 2
            else:
                H[i][j] = left; P[i][j] = 3

    # performing the backtracking based on the computed scoring matrix H
    j_star = max(range(0, m + 1), key=lambda jj: H[n][jj])
    best = H[n][j_star]
    i, j = n, j_star
    a1, a2 = [], []
    while i > 0 and j > 0 and P[i][j] != 0:
        if P[i][j] == 1:
            a1.append(A[i - 1]); a2.append(B[j - 1])
            i -= 1
            j -= 1
        elif P[i][j] == 2:
            a1.append(A[i - 1]); a2.append("-");
            i -= 1
        else:
            a1.append("-");      a2.append(B[j - 1])
            j -= 1
    
    a1.reverse()
    a2.reverse()

    # returning the based alignment, the alignment score and the indexes of the alignment
    return AlignmentRecord(
        score=best, i_start=i+1, i_end=n, j_start=1, j_end=j_star,
        aln1="".join(a1),
        aln2="".join(a2)
    )

# Storing the contig sequence and the sequences that took part in order
@dataclass
class Contig:
    seq: str
    ids: []

# Exact K-mer comparison by the suffix (hashing, so O(1) time)
def exact_kmer_right(a, b, k):
    return len(a) >= k and len(b) >= k and a[-k:] == b[:k]

# Exact K-mer comparison by the prefix (hashing, so O(1) time)
def exact_kmer_left(a, b, k):
    return len(a) >= k and len(b) >= k and b[-k:] == a[:k]

# Deterministic neighbor selection
# Scans for the best k-mer matching pair, which will be used for the DP extension
def find_best_extension_for(contig, others):
    candidates = []
    # RIGHT: contig.suffix -> other.prefix
    for j, o in enumerate(others):
        if exact_kmer_right(contig.seq, o.seq, KMER):
            rec = overlap_alignment_suffix_prefix(contig.seq, o.seq)
            if rec.j_end >= KMER:
                key = (rec.score, 1, min(o.ids))
                candidates.append((j, "RIGHT", rec, key))

    # LEFT: other.suffix -> contig.prefix
    for j, o in enumerate(others):
        if exact_kmer_left(contig.seq, o.seq, KMER):
            rec = overlap_alignment_suffix_prefix(o.seq, contig.seq)
            if rec.j_end >= KMER:
                key = (rec.score, 0, min(o.ids))
                candidates.append((j, "LEFT", rec, key))

    if not candidates:
        return None
    
    best_j, best_dir, best_rec, _ = max(candidates, key=lambda x: x[3])

    return (best_j, best_dir, best_rec)

# Merge two sequences with prefix/suffix match (extension)
def merge_pair(contig, partner, direction, rec):
    if direction == "RIGHT":
        new_seq = contig.seq + partner.seq[KMER:]
        new_ids = contig.ids + partner.ids
    else:  # LEFT
        new_seq = partner.seq + contig.seq[KMER:]
        new_ids = partner.ids + contig.ids

    return Contig(seq=new_seq, ids=new_ids)

# Deterministic seed selection for the k-mer match and extension by DP. Done in lexicographical order.
def _seed_index(cands):
    return min(range(len(cands)), key=lambda i: (min(cands[i].ids), len(cands[i].ids), len(cands[i].seq)))

# Main assembly process
def classical_assembly(reads):
    candidates = [Contig(seq=s, ids=[rid]) for rid, s in reads]
    contigs = []

    while candidates:
        idx = _seed_index(candidates)
        contig = candidates.pop(idx)
        if VERBOSE:
            print(f"[seed] {min(contig.ids)} (len={len(contig.seq)})")

        extended = True
        while extended and candidates:
            extended = False
            best = find_best_extension_for(contig, candidates)
            if best is None:
                break
            partner_idx, direction, rec = best
            partner = candidates.pop(partner_idx)
            if VERBOSE:
                left_id  = (contig.ids[-1] if direction == "RIGHT" else partner.ids[-1])
                right_id = (partner.ids[0]   if direction == "RIGHT" else contig.ids[0])
                print(f"  [{direction}] {left_id} -> {right_id}  (score={rec.score}, overlap>=K={KMER})")
            contig = merge_pair(contig, partner, direction, rec)
            extended = True

        if VERBOSE:
            print(f"[assembled contig] len={len(contig.seq)}  order={'->'.join(contig.ids)}")
        contigs.append(contig)

    return contigs

# setting the contid id and information of its assembly (includes order of assembled seqs)
def contigs_to_fasta_records(contigs):
    records = []
    for i, c in enumerate(contigs, start=1):
        header = f"contig_{i} nreads={len(c.ids)} len={len(c.seq)} order={'->'.join(c.ids)}"
        records.append((header, c.seq))
    return records



# main execution
def main():
    print(f"=== Running the Classical DP based Assembler ===")
    print(f"[project] {PROJECT_DIR}")
    print(f"[data]    {DATA_DIR}")
    print(f"[input]   {IN_FASTA}")
    print(f"[res]     {RES_DIR}")

    # read input FASTA file
    reads = read_fasta(IN_FASTA)
    
    # perform the assembly
    print("\n--- Assembly START ---")
    contigs = classical_assembly(reads)
    print("\n--- Assembly END ---")

    # store assembed result in FASTA format
    fasta_records = contigs_to_fasta_records(contigs)

    # printing assembly results to stdout
    print("\n--- Assembled Contig ---")
    print_fasta(fasta_records, line_wrap=LINE_WRAP)

    # write assembled result to file
    out_fasta = RES_DIR / "sarscov_S_10mer_seqs_shuffled_classical_assembled.fasta"
    write_fasta(fasta_records, out_fasta, line_wrap=LINE_WRAP)

    # printing summary and assembly order
    print("\n--- Assembly summary and order ---")
    print(f"[contigs] assembled contig written to {out_fasta}")
    for i, c in enumerate(contigs, start=1):
        print(f"  - contig_{i}: len={len(c.seq)} reads={len(c.ids)}")
        print(f"  - order={'->'.join(c.ids)}")

if __name__ == "__main__":
    main()