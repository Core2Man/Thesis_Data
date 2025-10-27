# ProMap

**ProMap** is a lightweight Python tool for mapping bacterial protein identifiers (RefSeq, UniProt, TRO) and retrieving predicted subcellular localizations from DeepLocPro results.

---

## 🧬 Features

- Maps protein sequences across:
  - NCBI RefSeq
  - UniProt
  - TRO (custom/internal IDs)
- Falls back to local BLAST if no exact match is found
- Integrates localization predictions from DeepLocPro
- Works entirely offline with local databases
- Interactive prompt for entering protein IDs

---

## 📁 Files & Requirements

Place these files in the same directory:

- `promap.py` — main script
- `DeepLocDA.csv` — CSV file with `ACC` and `Localization` columns
- `*.fasta` files — RefSeq, UniProt, and TRO protein FASTA files
- BLAST databases built from those FASTA files

---

## 🔧 Setup

1. Create a virtual environment *(optional but recommended)*:
   ```bash
   python -m venv venv
   source venv/bin/activate  # or venv\Scripts\activate on Windows
