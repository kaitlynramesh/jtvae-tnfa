"""
TNF-alpha binding oracle based on Tanimoto similarity to known binders,
weighted by binding affinity. Higher score = more likely to bind TNF-alpha.

Score ∈ [0, 1]:
  - Tanimoto similarity to nearest known binder (weighted by affinity)
  - Boosted by drug-likeness (QED)
"""
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, QED
import os

DATA_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "data", "tnfa_jtvae", "binders_with_affinity.csv"
)

def _load_reference(path, max_affinity_nM=1000.0):
    """Load known binders with affinity <= threshold (potent binders only)."""
    df = pd.read_csv(path)
    df = df[df["affinity_nM"] <= max_affinity_nM].copy()
    fps, weights = [], []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        fps.append(fp)
        # Convert affinity to weight: lower nM → higher weight (pIC50-like)
        # pKi = -log10(Ki_M) = -log10(Ki_nM * 1e-9) = 9 - log10(Ki_nM)
        pki = 9.0 - np.log10(max(row["affinity_nM"], 1e-3))
        weights.append(pki)
    weights = np.array(weights)
    weights = (weights - weights.min()) / (weights.max() - weights.min() + 1e-8)
    return fps, weights


_REF_FPS, _REF_WEIGHTS = None, None

def _get_reference():
    global _REF_FPS, _REF_WEIGHTS
    if _REF_FPS is None:
        _REF_FPS, _REF_WEIGHTS = _load_reference(DATA_PATH)
    return _REF_FPS, _REF_WEIGHTS


NOVELTY_THRESHOLD = 0.4  # Tanimoto cutoff: max_sim < threshold → novel


def score_components(smiles: str) -> dict:
    """Return all score components for ablation analysis.

    Keys:
      weighted_sim  — affinity-weighted Tanimoto to top-5 reference binders (0-1)
      max_sim       — raw max Tanimoto to any reference binder (0-1)
      novelty       — 1 - max_sim; higher = more structurally distinct (0-1)
      is_novel      — 1 if max_sim < NOVELTY_THRESHOLD (0.4), else 0
      qed           — RDKit QED drug-likeness (0-1)
      final         — combined oracle score: 0.7*weighted_sim + 0.3*qed (0-1)
    """
    mol = Chem.MolFromSmiles(smiles)
    null = {"weighted_sim": 0.0, "max_sim": 0.0, "novelty": 1.0, "is_novel": 1, "qed": 0.0, "final": 0.0}
    if mol is None:
        return null
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    ref_fps, ref_weights = _get_reference()
    if not ref_fps:
        return null

    sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, ref_fps))
    top_k = 5
    top_idx = np.argsort(sims)[-top_k:]
    weighted_sim = float(np.sum(sims[top_idx] * ref_weights[top_idx]) / (np.sum(ref_weights[top_idx]) + 1e-8))
    max_sim = float(sims.max())
    novelty = float(1.0 - max_sim)
    is_novel = int(max_sim < NOVELTY_THRESHOLD)

    try:
        qed_score = float(QED.qed(mol))
    except Exception:
        qed_score = 0.5

    final = float(np.clip(0.7 * weighted_sim + 0.3 * qed_score, 0.0, 1.0))
    return {
        "weighted_sim": weighted_sim, "max_sim": max_sim,
        "novelty": novelty, "is_novel": is_novel,
        "qed": qed_score, "final": final,
    }


def score(smiles: str) -> float:
    """Score a SMILES for predicted TNF-alpha binding. Returns float in [0, 1]."""
    return score_components(smiles)["final"]


def batch_score(smiles_list):
    return [score(s) for s in smiles_list]


def batch_score_components(smiles_list):
    """Return a DataFrame with one row per molecule and columns for each component."""
    rows = [score_components(s) for s in smiles_list]
    return pd.DataFrame(rows)


if __name__ == "__main__":
    test_smiles = [
        "COc1ccc2c(c1)C(=O)N(C[C@@]1(C#Cc3cncc(-c4cn[nH]c4)c3)NC(=O)NC1=O)C2",     # top binder
        "c1ccccc1",                                                                # benzene (low)
        "CC(=O)Oc1ccccc1C(=O)O",                                                   # aspirin
    ]
    for smi in test_smiles:
        print(f"Score: {score(smi):.4f}  |  {smi[:50]}")
