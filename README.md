# 🧬 Quantum Assembler (Classical + Grover)

Two proof-of-concept (PoC) assemblers for **k-mer overlap assembly**:

- **classical_assembler.py** — classical dynamic programming (DP)-based assembler (no quantum)
- **grover_assembler.py** — Grover-accelerated DP assembler implemented with the **Qiskit Aer** simulator

This project explores how **quantum search (Grover’s algorithm)** can be adapted to accelerate parts of the genome assembly process.

---

## 📁 Project Layout

```text
package/
├─ bin/
│  ├─ classical_assembler.py
│  └─ grover_assembler.py
├─ data/
│  ├─ sarscov_N_4mer_shuffled.fasta      # k-merized sequences of SARS-CoV-2 N gene (K=4)
│  ├─ sarscov_N_10mer_shuffled.fasta     # k-merized sequences of SARS-CoV-2 N gene (K=10)
│  └─ sarscov_N.fasta                    # full-length SARS-CoV-2 N gene
├─ res/
│  ├─ sarscov_S_4mer_seqs_shuffled_grover_assembled.fasta       # Grover-DP assembler output
│  └─ sarscov_S_10mer_seqs_shuffled_classical_assembled.fasta   # Classical assembler output
├─ setup.py      # installation script for dependencies (Qiskit, Aer, etc.)
└─ README.md


⸻

⚙️ Installation

# Clone this repository
git clone https://github.com/<your-username>/<your-repo-name>.git
cd <your-repo-name>

# Create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Upgrade pip and install dependencies via setup.py
pip install --upgrade pip
pip install .

💡 The setup.py script automatically installs qiskit and qiskit-aer, required for Grover’s quantum simulation.

⸻

🚀 Usage

Run the Classical Assembler

cd package/bin
python classical_assembler.py

Run the Grover-DP Assembler

cd package/bin
python grover_assembler.py

Both scripts assemble shuffled k-mers from the provided FASTA files and write the reconstructed sequence to the res/ directory.

⸻

🧠 Concepts
	•	k-mer assembly: reconstructing a full DNA sequence from overlapping k-length substrings.
	•	Grover’s algorithm: a quantum search algorithm that finds a marked item in an unsorted database in O(√N) time.
	•	Grover-DP hybrid: combines Grover’s amplitude amplification with classical dynamic programming (DP) scoring to search optimal overlaps between k-mers.

⸻

📊 Example Data

The project includes small test datasets derived from the SARS-CoV-2 N gene, preprocessed into shuffled k-mers:
	•	data/sarscov_N_4mer_shuffled.fasta
	•	data/sarscov_N_10mer_shuffled.fasta

You can substitute your own sequences by placing FASTA files in the data/ directory.

⸻

🧩 Output

Results are written under res/, for example:

res/
├─ sarscov_S_4mer_seqs_shuffled_grover_assembled.fasta       # Quantum version
└─ sarscov_S_10mer_seqs_shuffled_classical_assembled.fasta   # Classical version


⸻

🪐 Requirements
	•	Python ≥ 3.10
	•	Qiskit ≥ 1.0
	•	Qiskit Aer

Install them manually if needed:

pip install qiskit qiskit-aer


⸻

🧾 License

This project is licensed under the MIT License — free to use, modify, and distribute with attribution.

Copyright (c) 2025 Inuk Jung and Yongtae Kim
Kyungpook National University

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
[...]


⸻

👥 Authors

Inuk Jung and Yongtae Kim
Kyungpook National University (KNU), 2025
🧩 Quantum-Assisted DNA Assembly Project

⸻

🌌 Citation

If you use this code in your work, please cite:

Jung I, Kim Y. Quantum-Assisted Genome Assembly using Grover’s Algorithm.
Kyungpook National University, 2025.


⸻

💬 Contact

For questions or collaborations:
	•	📧 inukjung@knu.ac.kr
	•	🌐 https://github.com/

⸻

“Where quantum superposition meets biological reconstruction.” — Quantum Assembler

---
