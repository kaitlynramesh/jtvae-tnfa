"""
Run JT-VAE Bayesian Optimization to generate candidate TNF-alpha binders.

Uses the pretrained ZINC JT-VAE model as the molecular generator and a
TNF-alpha Tanimoto-similarity oracle as the scoring function.

Usage:
    conda activate jtvae-tnfa
    python run_jtvae_tnfa.py [--n_iters N] [--bo_batch N] [--n_seed N] [--seed N]
"""

import os
import sys
import argparse
import random
import numpy as np
import torch
import pandas as pd
from datetime import datetime

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
JT_VAE_DIR  = os.path.join(SCRIPT_DIR, "jt_vae")

sys.path.insert(0, JT_VAE_DIR)
sys.path.insert(0, SCRIPT_DIR)

import rdkit
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from fast_jtnn import JTNNVAE, Vocab
from tnfa_oracle import score as tnfa_score

from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_model
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import UpperConfidenceBound
from botorch.optim import optimize_acqf


# ── Paths ──────────────────────────────────────────────────────────────────────
VOCAB_PATH  = os.path.join(JT_VAE_DIR, "data", "zinc", "vocab.txt")
MODEL_PATH  = os.path.join(JT_VAE_DIR, "fast_molvae", "vae_model", "model.iter-25000")
SEED_SMILES = os.path.join(SCRIPT_DIR, "data", "tnfa_jtvae", "top100_binders.txt")
OUT_DIR     = os.path.join(SCRIPT_DIR, "results", "jtvae_tnfa")


def load_vae(vocab_path, model_path, device):
    vocab = Vocab([x.strip() for x in open(vocab_path)])
    model = JTNNVAE(vocab, hidden_size=450, latent_size=56, depthT=20, depthG=3)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model = model.to(device)
    model.eval()
    return model


def encode(model, smiles_list, device):
    """Encode a list of SMILES to latent vectors; skip failures."""
    vecs, valid = [], []
    for smi in smiles_list:
        try:
            z = model.encode_latent_mean([smi])
            vecs.append(z.detach().cpu())
            valid.append(smi)
        except Exception:
            pass
    if not vecs:
        return None, []
    return torch.cat(vecs, dim=0), valid


def decode_candidates(model, z_batch, n_decode=3):
    """Decode latent vectors to SMILES; try n_decode times for diversity."""
    candidates = set()
    device = next(model.parameters()).device
    for z in z_batch:
        z = z.view(1, -1).to(device)
        tree_vec, mol_vec = torch.hsplit(z, 2)
        for prob in [True, False] * ((n_decode + 1) // 2):
            try:
                smi = model.decode(tree_vec, mol_vec, prob_decode=prob)
                if smi:
                    candidates.add(smi)
            except Exception:
                pass
    return list(candidates)


def run_bo(model, seed_smiles, n_iters, bo_batch, device):
    """Bayesian Optimisation loop in the JT-VAE latent space."""
    print(f"Encoding {len(seed_smiles)} seed molecules...")
    train_X, valid_seed = encode(model, seed_smiles, device)
    if train_X is None or len(valid_seed) == 0:
        raise RuntimeError("No seed SMILES could be encoded – check vocab/model.")

    train_Y = torch.tensor([tnfa_score(s) for s in valid_seed]).view(-1, 1).float()
    print(f"  Encoded {len(valid_seed)} seeds. Score range: "
          f"{train_Y.min():.3f} – {train_Y.max():.3f}")

    all_results = []
    for s, y in zip(valid_seed, train_Y.squeeze().tolist()):
        all_results.append({"smiles": s, "score": y, "iteration": 0})

    bounds = torch.stack([train_X.min(0).values, train_X.max(0).values])

    for iteration in range(1, n_iters + 1):
        print(f"\n── Iteration {iteration}/{n_iters} ──")

        # 1. Fit GP
        gp = SingleTaskGP(train_X, train_Y)
        mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
        fit_gpytorch_model(mll)

        # 2. Optimise acquisition (UCB) — q=1 per botorch UCB constraint; loop bo_batch times
        ucb = UpperConfidenceBound(gp, beta=0.1)
        candidates_z = []
        for _ in range(bo_batch):
            z_cand, _ = optimize_acqf(
                ucb,
                bounds=bounds,
                q=1,
                num_restarts=5,
                raw_samples=20,
            )
            candidates_z.append(z_cand.squeeze(0))
        candidates_z = torch.stack(candidates_z)

        # 3. Decode candidates
        new_smiles = decode_candidates(model, candidates_z, n_decode=3)
        print(f"  Decoded {len(new_smiles)} candidate molecules")

        if not new_smiles:
            print("  No valid molecules decoded – skipping iteration.")
            continue

        # 4. Score and encode new molecules
        new_X_list, new_Y_list, scored = [], [], []
        for smi in new_smiles:
            try:
                z = model.encode_latent_mean([smi]).detach().cpu()
                s = tnfa_score(smi)
                new_X_list.append(z)
                new_Y_list.append(s)
                scored.append({"smiles": smi, "score": s, "iteration": iteration})
            except Exception:
                pass

        if not new_X_list:
            continue

        new_X = torch.cat(new_X_list, dim=0)
        new_Y = torch.tensor(new_Y_list).view(-1, 1).float()

        # 5. Update training data
        train_X = torch.cat([train_X, new_X], dim=0)
        train_Y = torch.cat([train_Y, new_Y], dim=0)
        all_results.extend(scored)

        best_iter = max(new_Y_list)
        best_all  = train_Y.max().item()
        print(f"  Best this iter: {best_iter:.4f} | Best overall: {best_all:.4f}")

    return all_results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_iters",  type=int, default=20,
                        help="Number of BO iterations (default: 20)")
    parser.add_argument("--bo_batch", type=int, default=10,
                        help="Candidates proposed per BO step (default: 10)")
    parser.add_argument("--n_seed",   type=int, default=50,
                        help="Seed molecules from top binders (default: 50)")
    parser.add_argument("--seed",     type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    os.makedirs(OUT_DIR, exist_ok=True)

    print("Loading JT-VAE model...")
    model = load_vae(VOCAB_PATH, MODEL_PATH, device)
    print("Model loaded.")

    seed_smiles = [l.strip() for l in open(SEED_SMILES) if l.strip()]
    seed_smiles = seed_smiles[:args.n_seed]
    print(f"Using {len(seed_smiles)} seed molecules.")

    results = run_bo(model, seed_smiles, args.n_iters, args.bo_batch, device)

    df = pd.DataFrame(results).sort_values("score", ascending=False).drop_duplicates("smiles")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(OUT_DIR, f"candidates_{timestamp}.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved {len(df)} candidate molecules → {out_path}")

    print("\nTop 10 generated candidates:")
    print(df.head(10)[["smiles", "score", "iteration"]].to_string(index=False))


if __name__ == "__main__":
    main()
