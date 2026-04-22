# JT-VAE for TNF-alpha Binder Generation

CBB 5801 Proteomics / Protein Modeling Group Project

We use the [mol_opt](https://github.com/wenhao-gao/mol_opt) implementation of Junction Tree Variational Autoencoder (JT-VAE) to generate small-molecule binder candidates for TNF-alpha. This version of JT-VAE is pretrained on the ZINC-250K database, so we refine this model on known TNF-alpha binders from BindingDB.

## Set up environment
```bash
git clone https://github.com/kaitlynramesh/jtvae-tnfa.git
cd jtvae-tnfa
bash setup.sh # set up conda environment using yml file 
```

The earlier version of JT-VAE in mol_opt uses GPU, but this implementation is only compatible with CPU. Switching to GPU can be done in the future. 

## Overview
```
conda activate jtvae-tnfa
python prepare_data.py 
python run_jtvae_tnfa.py --n_iters 20 --bo_batch 10 --n_seed 50 --seed 42
# n_seed = number of "seed" molecules taken from known binders
# seed = seed for reproducibility
```
Results stored in `results/jtvae_tnfa_candidates_*.csv`

## Scoring function 
See `tnfa_oracle.py`
- Similarity score: Weighted sum of Tanimoto similarity to known TNF-alpha binders. The "weight" is set by the K_i or IC50 value corresponding to the nearest binder. This weighted sum based on the top five most similar binders are used to score each candidate molecule.
- QED (drug-likeness score) that accounts for chemical and structural properties of the moelcule.
Final score = 0.7 * similarity + 0.3 * QED

We can check the individual score components by running
```
from tnfa_oracle import score_components
score_components(<SMILES string>)
```
## Thresholding
The candidate sequences of interest for analysis with molecular docking are filtered based on novelty (Tanimoto similarity <= 0.4)

## References
- Claude Code 3.5 Sonnet
- **JT-VAE**: Jin et al., "Junction Tree Variational Autoencoder for Molecular Graph Generation," ICML 2018. [arXiv:1802.04364](https://arxiv.org/abs/1802.04364)
- **mol_opt benchmark**: Gao et al., "Sample Efficiency Matters: A Benchmark for Practical Molecular Optimization," NeurIPS 2022. [GitHub](https://github.com/wenhao-gao/mol_opt)
- **AutoDock Vina (original)**: Trott & Olson, "AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading," J. Comput. Chem. 2010, 31(2):455-461. DOI
- **AutoDock Vina**: Eberhardt et al., "AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings," J. Chem. Inf. Model. 2021, 61(8):3891-3898. DOI
  
