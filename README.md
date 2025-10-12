# Quantum Assembler (Classical + Grover)

Two PoC assemblers for k-mer overlap assembly:

- classical_assembler.py — the classical DP based assembler (no quantum)
- grover_assembler.py — the Grover-DP assembler (using Qiskit Aer simulator)


---

## Project layout
package/
├─ bin/
│  ├─ classical_assembler.py
│  └─ grover_assembler.py
├─ data/
│  ├─ sarscov_N_4mer_shuffled.fasta     # k-merized sequences of the SARS-CoV-2 N gene (K=4)
│  └─ sarscov_N_10mer_shuffled.fasta    # k-merized sequences of the SARS-CoV-2 N gene (K=10)
│  └─ sarscov_N.fasta                   # full length sequence of the SARS-CoV-2 N gene
├─ res/
│  └─ sarscov_S_4mer_seqs_shuffled_grover_assembled.fasta       # output of Grover-DP assembler
│  └─ sarscov_S_10mer_seqs_shuffled_classical_assembled.fasta   # output of classical assembler
├─ setup.py     # script to install qskit packages
└─ README.md

## Install

```bash
# from the project root (where setup.py lives)
python3 -m venv .venv
source .venv/bin/activate

# install dependencies via setup.py
pip install --upgrade pip
pip install .

# run the classical assembler
cd bin
python classical_assembler.py

# run the Grover-DP assembler
cd bin
python grover_assembler.py