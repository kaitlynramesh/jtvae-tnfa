# JT-VAE TNF-alpha Binder Generation

Uses the [mol_opt](https://github.com/wenhao-gao/mol_opt) framework's Junction Tree Variational Autoencoder (JT-VAE) with Bayesian Optimization to generate candidate small-molecule binders for **TNF-alpha** (Tumor Necrosis Factor alpha), starting from known binders in BindingDB.

---

## Overview

The pipeline has three stages:

```
BindingDB TSV  ──►  prepare_data.py  ──►  known binder SMILES + affinities
                                                       │
                                                       ▼
                                             run_jtvae_tnfa.py
                                        (JT-VAE encoder → GP-BO → decoder)
                                                       │
                                                       ▼
                                          results/jtvae_tnfa/candidates_*.csv
```

**JT-VAE** encodes molecules into a continuous 112-dimensional latent space. **Bayesian Optimization** (Gaussian Process + Upper Confidence Bound acquisition) navigates that space to find new molecules predicted to score highly. Candidate latent vectors are decoded back to SMILES.

The **oracle** (scoring function) combines:
- Weighted Tanimoto similarity to known TNF-alpha binders (70%) — potent binders (low Ki/IC50) are weighted more heavily
- QED drug-likeness score (30%)

---

## Repository Layout

```
jtvae-tnfa/
├── jt_vae/                                  # JT-VAE code (from wenhao-gao/mol_opt, patched for CPU)
│   ├── fast_jtnn/                           # model architecture
│   ├── fast_molvae/
│   │   └── vae_model/model.iter-25000       # pretrained ZINC VAE weights (21 MB)
│   └── data/zinc/vocab.txt                  # ZINC vocabulary
├── data/
│   ├── TNFA_CBB_FINAL_PROJECT.tsv           # raw BindingDB export
│   └── tnfa_jtvae/
│       ├── all_binders.txt                  # 1,354 canonicalized unique SMILES
│       ├── top100_binders.txt               # top-100 by Ki/IC50 (BO seed)
│       └── binders_with_affinity.csv        # SMILES + best affinity_nM
├── results/
│   └── jtvae_tnfa/
│       └── candidates_<timestamp>.csv       # generated candidates (output)
├── prepare_data.py                          # Step 1: extract/clean SMILES from TSV
├── tnfa_oracle.py                           # Step 2: scoring function
├── run_jtvae_tnfa.py                        # Step 3: BO optimization runner
├── scoring_eval.ipynb                       # analysis notebook
├── environment.yml                          # reproducible conda environment
└── setup.sh                                 # one-command environment setup + verification
```

---

## Environment Setup

### Prerequisites

- [Anaconda or Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Git
- macOS or Linux (CPU-only; for CUDA see note below)

### Quick setup (recommended)

```bash
git clone https://github.com/<your-username>/jtvae-tnfa.git
cd jtvae-tnfa
bash setup.sh
```

`setup.sh` creates the `jtvae-tnfa` conda environment from `environment.yml`, then loads the pretrained model and decodes a test molecule to confirm everything works.

### Manual setup (if you prefer step-by-step)

```bash
conda env create -f environment.yml
conda activate jtvae-tnfa
```

> **CUDA:** The default `environment.yml` installs the CPU-only PyTorch wheel. For GPU, edit the `torch` line to your CUDA variant from [pytorch.org](https://pytorch.org/get-started/locally/).

### CPU compatibility patches (already applied)

The upstream JT-VAE code assumes a GPU. `jt_vae/fast_jtnn/` in this repo has two patches applied:

- **`nnutils.py` — `create_var`**: falls back to CPU tensors when `torch.cuda.is_available()` is False
- **`jtnn_vae.py` — `sample_prior` / `dfs_assemble`**: replaces hardcoded `.cuda()` with `.to(device)` calls

If you ever re-pull the upstream JT-VAE code, these patches need to be re-applied.

---

## Running the Pipeline

### Step 1 — Prepare data (run once)

Extracts and canonicalizes SMILES from the BindingDB TSV, filters to entries with binding data, deduplicates, and saves output files to `data/tnfa_jtvae/`.

```bash
conda activate molopt
cd /path/to/0416_MolOpt
python prepare_data.py
```

Expected output:
```
Total valid binders: 1354
Affinity range: 0.030 – 20600000.0 nM
Saved 1354 SMILES → data/tnfa_jtvae/all_binders.txt
Saved affinity data → data/tnfa_jtvae/binders_with_affinity.csv
Saved top-100 binders → data/tnfa_jtvae/top100_binders.txt
```

### Step 2 — Verify the oracle (optional)

```bash
python tnfa_oracle.py
```

Expected output:
```
Score: 0.7126  |  COc1ccc2c(c1)C(=O)N(...)   # top known binder — high score
Score: 0.1774  |  c1ccccc1                     # benzene — low score
Score: 0.3098  |  CC(=O)Oc1ccccc1C(=O)O       # aspirin — low score
```

To inspect all score components (similarity, novelty, QED) for a molecule:
```python
from tnfa_oracle import score_components
score_components("c1ccccc1")
# {'weighted_sim': 0.0, 'max_sim': 0.07, 'novelty': 0.93, 'is_novel': 1, 'qed': 0.36, 'final': 0.14}
```

### Step 3 — Run JT-VAE Bayesian Optimization

```bash
python run_jtvae_tnfa.py --n_iters 20 --bo_batch 10 --n_seed 50
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--n_iters` | 20 | Number of BO iterations |
| `--bo_batch` | 10 | New candidates proposed per iteration |
| `--n_seed` | 50 | Seed molecules from top known binders |
| `--seed` | 42 | Random seed for reproducibility |

Output is saved to `results/jtvae_tnfa/candidates_<timestamp>.csv` with columns:

| Column | Description |
|--------|-------------|
| `smiles` | Canonical SMILES of generated molecule |
| `score` | Oracle score ∈ [0, 1] (higher = better predicted binder) |
| `iteration` | BO iteration that produced it (0 = seed molecules) |

---

## How JT-VAE + Bayesian Optimization Works

```
Known binders (SMILES)
        │
        ▼
  JT-VAE Encoder ──► Latent vectors z ∈ ℝ¹¹² (56 tree + 56 graph)
        │
        ▼
  Score each z with oracle
        │
        ▼
  Fit Gaussian Process on (z, score) pairs
        │
        ▼
  Maximize UCB acquisition to find promising z*
        │
        ▼
  JT-VAE Decoder (z* → SMILES)
        │
        ▼
  Score new molecule, add to GP training set
        │
        └──► repeat for n_iters
```

The JT-VAE represents molecules as junction trees (rings and chains as nodes) assembled from a vocabulary of molecular fragments. This guarantees all decoded molecules are **chemically valid**.

---

## Oracle Details

The TNF-alpha oracle in `tnfa_oracle.py` does not use docking — it uses a ligand-based proxy:

1. **Morgan fingerprints** (radius=2, 2048 bits) computed for the query molecule
2. **Tanimoto similarity** computed against all known potent binders (Ki or IC50 ≤ 1000 nM)
3. Each reference binder is **weighted by pKi** (`9 − log10(Ki_nM)`) so tight binders contribute more
4. **Top-5 nearest neighbour** weighted similarity is taken as the binding signal
5. **QED** (Quantitative Estimate of Drug-likeness) adds a 30% drug-likeness bonus

Final score = `0.7 × weighted_similarity + 0.3 × QED`

> This is a heuristic oracle. For downstream validation, generated candidates should be assessed with molecular docking (e.g., AutoDock Vina against PDB structure 2AZ5) or experimental assays.

### Score components API

`score_components(smiles)` returns a dict with all intermediate values, useful for ablation analysis:

| Key | Description |
|-----|-------------|
| `weighted_sim` | Affinity-weighted Tanimoto to top-5 reference binders (0–1) |
| `max_sim` | Raw maximum Tanimoto to any reference binder (0–1) |
| `novelty` | `1 − max_sim`; higher = more structurally distinct from known binders (0–1) |
| `is_novel` | `1` if `max_sim < 0.4` (standard Tanimoto novelty cutoff), else `0` |
| `qed` | RDKit QED drug-likeness (0–1) |
| `final` | Combined oracle score: `0.7 × weighted_sim + 0.3 × qed` (0–1) |

`batch_score_components(smiles_list)` returns a pandas DataFrame with the same columns.

The novelty threshold (`NOVELTY_THRESHOLD = 0.4`) is exported from the module so the evaluation notebook and oracle stay in sync.

---

## Data Source

| Field | Value |
|-------|-------|
| Database | [BindingDB](https://www.bindingdb.org/) |
| Target | Tumor Necrosis Factor alpha (TNF-α), UniProt P01375 |
| File | `TNFA_CBB_FINAL_PROJECT.tsv` |
| Entries | 3,163 rows → 1,354 unique valid binders with affinity data |
| Affinity range | 0.030 nM – 20,600,000 nM (Ki or IC50) |

---

## Key References

- **JT-VAE**: Jin et al., "Junction Tree Variational Autoencoder for Molecular Graph Generation," ICML 2018. [arXiv:1802.04364](https://arxiv.org/abs/1802.04364)
- **mol_opt benchmark**: Gao et al., "Sample Efficiency Matters: A Benchmark for Practical Molecular Optimization," NeurIPS 2022. [GitHub](https://github.com/wenhao-gao/mol_opt)
- **BoTorch**: Balandat et al., "BoTorch: A Framework for Efficient Monte-Carlo Bayesian Optimization," NeurIPS 2020.
