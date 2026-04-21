#!/usr/bin/env bash
# One-time setup: create conda environment and verify the installation.
set -e

echo "=== Creating conda environment from environment.yml ==="
conda env create -f environment.yml
conda activate jtvae-tnfa

echo ""
echo "=== Verifying installation ==="
python - <<'EOF'
import torch
from rdkit import Chem
from fast_jtnn import JTNNVAE, Vocab
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "jt_vae"))

vocab_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "jt_vae", "data", "zinc", "vocab.txt")
model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "jt_vae", "fast_molvae", "vae_model", "model.iter-25000")

vocab = Vocab([x.strip() for x in open(vocab_path)])
model = JTNNVAE(vocab, hidden_size=450, latent_size=56, depthT=20, depthG=3)
model.load_state_dict(torch.load(model_path, map_location="cpu"))
model.eval()
print("JT-VAE model loaded successfully.")

smi = model.sample_prior()
print(f"Sample decoded SMILES: {smi}")

mol = Chem.MolFromSmiles(smi)
print(f"RDKit parse: {'valid' if mol else 'INVALID'}")
print("Setup complete — ready to run run_jtvae_tnfa.py")
EOF
