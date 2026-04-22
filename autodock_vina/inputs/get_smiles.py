import os
import glob
import pandas as pd
import subprocess
from rdkit import Chem

# -------------------------
# PATHS
# -------------------------
INPUT_DIR = "/Users/stuytschaevers/Documents/Yale/spring2026/comp_proteomics/final_project/final_docking_files/test_ligands_docking_output"
OUTPUT_CSV = "/Users/stuytschaevers/Documents/Yale/spring2026/comp_proteomics/final_project/final_docking_files/smiles_output/pdbqt_smiles.csv"

os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)

OBABEL = "obabel"


# -------------------------
# STEP 1: PDBQT → SMILES (OpenBabel ONLY)
# -------------------------
def pdbqt_to_smiles(pdbqt_path):
    cmd = [OBABEL, pdbqt_path, "-osmi"]

    p = subprocess.run(cmd, capture_output=True, text=True)

    if p.returncode != 0:
        return None, p.stderr.strip()

    lines = p.stdout.strip().splitlines()
    if not lines:
        return None, "empty output"

    smiles = lines[0].split()[0]

    return smiles, None


# -------------------------
# STEP 2: CLEAN + VALIDATE (RDKit ONLY)
# -------------------------
def clean_smiles(smiles):
    if not smiles:
        return None

    # parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # sanitize + canonicalize
    try:
        Chem.SanitizeMol(mol)
    except:
        return None

    return Chem.MolToSmiles(mol, canonical=True)


# -------------------------
# STEP 3: PIPELINE
# -------------------------
def process_file(pdbqt_file):
    raw_smiles, err = pdbqt_to_smiles(pdbqt_file)

    if not raw_smiles:
        return None, f"obabel error: {err}"

    clean = clean_smiles(raw_smiles)

    if not clean:
        return None, "RDKit failed to parse"

    return clean, None


# -------------------------
# MAIN
# -------------------------
def main():
    pdbqt_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.pdbqt")))

    rows = []

    for f in pdbqt_files:
        smiles, err = process_file(f)

        rows.append({
            "file_name": os.path.basename(f),
            "file_path": f,
            "smiles": smiles,
            "status": "OK" if smiles else "FAIL",
            "error": err
        })

        print(os.path.basename(f), "->", smiles if smiles else f"[FAIL] {err}")

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_CSV, index=False)

    print("\nWrote:", OUTPUT_CSV)


if __name__ == "__main__":
    main()
