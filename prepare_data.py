"""
Extract TNF-alpha binder SMILES from BindingDB TSV and prepare input files for jt-vae.
"""
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TSV_PATH = os.path.join(DATA_DIR, "TNFA_CBB_FINAL_PROJECT.tsv")
OUT_DIR = os.path.join(DATA_DIR, "tnfa_jtvae")

def canonicalize(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    df = pd.read_csv(TSV_PATH, sep="\t", header=0, low_memory=False)

    # Keep rows with SMILES and at least one binding measurement
    has_data = df["Ki (nM)"].notna() | df["IC50 (nM)"].notna()
    df = df[has_data & df["Ligand SMILES"].notna()].copy()

    # Canonical SMILES
    df["canonical_smiles"] = df["Ligand SMILES"].apply(canonicalize)
    df = df.dropna(subset=["canonical_smiles"])

    # Best affinity: prefer Ki, fall back to IC50 (both in nM)
    df["affinity_nM"] = df["Ki (nM)"].combine_first(df["IC50 (nM)"])
    df["affinity_nM"] = pd.to_numeric(df["affinity_nM"], errors="coerce")
    df = df.dropna(subset=["affinity_nM"])

    # De-duplicate: keep best affinity per SMILES
    df = df.sort_values("affinity_nM").drop_duplicates("canonical_smiles", keep="first")

    print(f"Total valid binders: {len(df)}")
    print(f"Affinity range: {df['affinity_nM'].min():.3f} – {df['affinity_nM'].max():.1f} nM")

    # Save all SMILES (one per line) for reference
    all_smi_path = os.path.join(OUT_DIR, "all_binders.txt")
    df["canonical_smiles"].to_csv(all_smi_path, index=False, header=False)
    print(f"Saved {len(df)} SMILES → {all_smi_path}")

    # Save with affinities for oracle training
    aff_path = os.path.join(OUT_DIR, "binders_with_affinity.csv")
    df[["canonical_smiles", "affinity_nM"]].to_csv(aff_path, index=False)
    print(f"Saved affinity data → {aff_path}")

    # Top-100 binders (lowest Ki) for seeding BO
    top_path = os.path.join(OUT_DIR, "top100_binders.txt")
    df.head(100)["canonical_smiles"].to_csv(top_path, index=False, header=False)
    print(f"Saved top-100 binders → {top_path}")

    print("\nTop 5 binders:")
    print(df[["canonical_smiles", "affinity_nM"]].head(5).to_string(index=False))

if __name__ == "__main__":
    main()
