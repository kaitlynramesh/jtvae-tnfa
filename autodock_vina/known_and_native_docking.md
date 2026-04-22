## README — Docking Top 5 Known TNFA Binders (by *K*ᵢ) + Native Ligand (AutoDock Vina)

### Overview
This workflow:
1. Loads a TSV of known TNFA ligands and their experimental *K*ᵢ values (nM).
2. Selects the **top 5** ligands with the **lowest *K*ᵢ** (strongest binders).
3. Converts each ligand SMILES to a **3D structure**, prepares it as **PDBQT** (Meeko), and docks it to a prepared receptor using **AutoDock Vina**.
4. Saves docking poses and a CSV of docking scores.
5. Separately docks a provided **native ligand** (from PDB 2AZ5) and writes its Vina score to disk.

---

### Inputs
Located under:
- `input_dir = .../final_docking_files/inputs`

Required files:
- `protein_chainA_B_fixed.pdbqt`  
  Prepared receptor in PDBQT format.
- `ad4_types.json`  
  Atom typing parameters for Meeko/AutoDock4-style types.
- `known_ligands/TNFA_CBB_FINAL_PROJECT.tsv`  
  Tab-separated table containing at least:
  - `Ki (nM)`
  - `Ligand SMILES`

---

### Outputs
Located under:
- `output_dir = .../final_docking_files/output`

Generated files/directories:
- `top5_known_database/ligand_top5_<i>.pdbqt`  
  Prepared ligand PDBQT for each of the top 5 ligands.
- `top5_known_database/docked_top5_<i>.pdbqt`  
  Docked poses for each ligand (top 5 poses written).
- `top5_known_docking_scores.csv`  
  Table containing (when successful): rank, smiles, vina_score, time_sec.
- `native/native_ligand.pdbqt`  
  Prepared native ligand.
- `native/native_docked_ligand.pdbqt`  
  Docked native ligand poses.
- `native/native_ligand_score.txt`  
  Text file containing the native ligand’s best-pose Vina score.

---

### Method Summary
#### 1) Ranking ligands by experimental affinity
- Reads the TSV (`sep="\t"`), sorts by `Ki (nM)` ascending.
- Takes the first 5 SMILES strings.

#### 2) SMILES → 3D → PDBQT (ligand preparation)
Function: `smiles_to_pdbqt(smiles, output_pdbqt, atom_params)`
- RDKit:
  - Parses SMILES
  - Adds hydrogens
  - Embeds 3D coordinates with ETKDG
  - Optimizes with UFF
- Meeko:
  - Assigns atom types/rotatable bonds
  - Writes PDBQT to disk

If RDKit cannot parse the SMILES, the ligand is skipped.

#### 3) Docking with AutoDock Vina
- Initializes Vina with `sf_name="vina"`.
- Sets receptor from `protein_pdbqt`.
- Computes maps once using:
  - `center=[-19.28, 74.56, 33.93]`
  - `box_size=[20, 20, 20]`
- For each ligand:
  - `exhaustiveness=16`
  - `n_poses=10` generated internally
  - Writes `n_poses=5` poses to file
  - Records the **top pose score** via `v.energies(n_poses=1)[0][0]`

#### 4) Native ligand docking
- Uses the hard-coded native SMILES:
  ```
  CN(CCN(C)Cc(coc1cc(C)c(C)cc12)c2=O)C[C@@H](c3ccccc43)CN4c5cc(C(F)(F)F)ccc5
  ```
- Prepares it similarly (RDKit + Meeko).
- Docks with the same grid and settings.
- Writes best score to `native_ligand_score.txt`.

---

### How to Run
1. Ensure the input paths in the script match your filesystem:
   - `input_dir`
   - `output_dir`
2. Confirm the receptor PDBQT and `ad4_types.json` exist in `input_dir`.
3. Run the notebook/script cells in order.

---

### Notes / Assumptions
- This pipeline assumes the receptor is already prepared and correctly protonated/typed for docking.
- The docking box center/size must correspond to the intended binding site.
- RDKit embedding/optimization can fail for some molecules; those are skipped (or may raise in the native-ligand section if invalid).
- Vina scores are in kcal/mol (more negative is better), and are not directly comparable to experimental *K*ᵢ without calibration.

---

### Key Dependencies
- Python: `numpy`, `pandas`
- Docking: `vina` (Python bindings)
- Chemistry: `rdkit`
- PDBQT prep: `meeko`

--- 

### Result Table (Top 5 Known Ligands)
The saved CSV `top5_known_docking_scores.csv` contains:
- `rank`: 1–5 in order of increasing experimental *K*ᵢ
- `smiles`: ligand SMILES
- `vina_score`: best docking score (top pose)
- `time_sec`: elapsed time per ligand docking run
