# 🧬 QuAssem — Quantum Assembler (Classical + Grover)

Two proof-of-concept assemblers for **k-mer overlap assembly**:

- **bin/classical_assembler.py** — classical DP-based assembler (no quantum)
- **bin/grover_assembler.py** — Grover-accelerated DP assembler using **Qiskit Aer**

This project explores how **Grover’s quantum search** can accelerate parts of genome assembly by amplifying overlap candidates that are then validated by classical dynamic programming (DP).

---

## 📁 Project Layout

```text
.
├─ bin/
│  ├─ classical_assembler.py
│  └─ grover_assembler.py
├─ data/
│  ├─ sarscov_N_4mer_shuffled.fasta      # k-merized sequences of SARS-CoV-2 N gene (K=4)
│  ├─ sarscov_N_10mer_shuffled.fasta     # k-merized sequences of SARS-CoV-2 N gene (K=10)
│  └─ sarscov_N.fasta                    # full-length SARS-CoV-2 N gene
├─ res/
│  ├─ sarscov_S_4mer_seqs_shuffled_grover_assembled.fasta       # Grover-DP output
│  └─ sarscov_S_10mer_seqs_shuffled_classical_assembled.fasta   # Classical output
├─ LICENSE
└─ README.md
```

> Tip: If you prefer a collapsible tree on GitHub, wrap the tree block with `<details><summary>Project tree</summary> ... </details>`.

---

## ⚙️ Installation

```bash
# clone
git clone https://github.com/cobi-git/QuAssem.git
cd QuAssem

# (optional) virtual env
python3 -m venv .venv
source .venv/bin/activate

# install deps
pip install --upgrade pip
pip install qiskit qiskit-aer
```

---

## 🚀 Usage

```bash
# classical assembler
python bin/classical_assembler.py

# Grover-DP assembler (Qiskit Aer)
python bin/grover_assembler.py
```

Outputs are written to `res/`.

---

## 🧠 Concepts (very brief)

- **k-mer assembly**: reconstruct a sequence from overlapping length-K substrings.  
- **Grover’s algorithm**: quantum O(√N) search for marked items; here, used to amplify valid overlap candidates.  
- **Hybrid Grover-DP**: quantum search proposes candidates; classical DP scoring validates/extends contigs.

---

## 🪐 Requirements

- Python ≥ 3.10
- Qiskit ≥ 1.0
- Qiskit Aer

Install manually if needed:
```bash
pip install qiskit qiskit-aer
```

---

## 🧾 License

**Apache License 2.0** — see `LICENSE` for details.

---

## 👥 Authors

Inuk Jung and Yongtae Kim  
Kyungpook National University (KNU), 2025
