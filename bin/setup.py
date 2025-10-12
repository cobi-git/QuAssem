#!/usr/bin/env python3

# Installs Python packages required by
#  - classical_assembler.py
#  - grover_assembler.py

import sys
import subprocess

PKGS = [
    "qiskit",
    "qiskit-aer",
]

def run(cmd):
    print(f"$ {' '.join(cmd)}", flush=True)
    subprocess.check_call(cmd)

def main():
    py = sys.executable

    try:
        run([py, "-m", "pip", "install", "--upgrade", "pip"])
    except subprocess.CalledProcessError:
        print("[warn] pip upgrade failed; continuing...", flush=True)

    try:
        run([py, "-m", "pip", "install", *PKGS])
        print("\n[ok] dependencies installed system/site-wide")
        return
    except subprocess.CalledProcessError:
        print("[warn] site install failed; retrying with --user ...", flush=True)

    run([py, "-m", "pip", "install", "--user", *PKGS])
    print("\n[ok] dependencies installed to user site-packages")

if __name__ == "__main__":
    try:
        main()
        print("\nAll set.")
    except subprocess.CalledProcessError as e:
        print(f"\n[error] command failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)