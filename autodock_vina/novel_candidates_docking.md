## README — Docking Top/Bottom Novel Candidates (by oracle_score) to TNFA (2AZ5 chains A/B)

### Overview
This script docks a subset of *novel candidate ligands* against TNFA using AutoDock Vina. It:
1. Prepares a receptor structure from a PDB file (chains A/B) using **PDBFixer** (adds missing atoms/residues and hydrogens).
2. Uses a pre-converted receptor **PDBQT** for docking.
3. Loads a CSV of candidate ligands with an `oracle_score`.
4. Selects the **top 10** and **bottom 10** candidates by `oracle_score` (deduplicated by SMILES).
5. Converts each ligand SMILES to 3D and prepares **ligand PDBQT** via **RDKit + Meeko**.
6. Docks each ligand with Vina, saves poses, and writes a CSV of docking scores.

---

### Inputs
Located under:
- `input_dir = .../final_docking_files/inputs`

Required files:
- `top_10_novel_candidates.csv`  
  Must contain columns:
  - `smiles`
  - `oracle_score`
- `2az5_chainA_B.pdb`  
  Receptor PDB containing only chains A and B (no waters/ligands/extra chains).
- `ad4_types.json`  
  AutoDock4 atom typing parameters used by Meeko.
- `protein_chainA_B_fixed.pdbqt`  
  Receptor in PDBQT format used by Vina (see “Receptor preparation + PDBQT conversion”).

Generated intermediate (written into `input_dir`):
- `protein_chainA_B_fixed.pdb`  
  Output of PDBFixer (receptor with missing atoms/residues and hydrogens added).

---

### Outputs
Located under:
- `output_dir = .../test_ligands_docking_output`

Generated files:
- `ligand_<idx>.pdbqt`  
  Prepared ligand file for each docked candidate (idx is the DataFrame index).
- `docked_<idx>.pdbqt`  
  Docked poses (top 5 poses written) for each candidate.
- `docking_scores.csv`  
  Summary table containing:
  - `ligand_id`
  - `ligand_file`
  - `smiles`
  - `vina_score` (best pose; kcal/mol; more negative is better)
  - `time_sec`

---

### Receptor preparation + PDBQT conversion
#### 1) PDBFixer receptor preparation (done in-script)
Function: `prepare_receptor(pdb_in, pdb_out)`
- Loads `2az5_chainA_B.pdb`
- Finds missing residues/atoms
- Adds missing atoms
- Adds hydrogens at **pH 7.4**
- Writes `protein_chainA_B_fixed.pdb`

#### 2) Convert receptor PDB → PDBQT (manual command-line step)
The script assumes you run (outside Python) something like:
- `obabel protein_fixed.pdb -O protein_fixed.pdbqt -xr`

Then docking uses:
- `protein_chainA_B_fixed.pdbqt`

---

### Ligand preparation (SMILES → 3D → PDBQT)
Function: `smiles_to_pdbqt(smiles, output_pdbqt, atom_params)`
- RDKit:
  - Parses SMILES
  - Adds hydrogens
  - 3D embedding (ETKDG)
  - UFF geometry optimization
- Meeko:
  - Assigns atom types/rotatable bonds using `ad4_types.json`
  - Writes ligand PDBQT

If SMILES parsing fails, that ligand is skipped.

---

### Docking configuration (AutoDock Vina)
- Scoring function: `sf_name="vina"`
- Grid maps computed once with:
  - `center = [-19.28, 74.56, 33.93]`
  - `box_size = [20, 20, 20]`
- For each ligand:
  - `exhaustiveness = 16`
  - `n_poses = 10` generated
  - `write_poses(..., n_poses=5)` saved
  - Best score read via `v.energies(n_poses=1)[0][0]`

---

### Candidate selection logic
1. Read `top_10_novel_candidates.csv` into a DataFrame.
2. Convert `oracle_score` to numeric; drop rows where it’s missing.
3. Sort by `oracle_score` descending (highest first).
4. Select:
   - `TOP_K = 10` highest-scoring
   - `TOP_K = 10` lowest-scoring
5. Concatenate both sets and `drop_duplicates(subset="smiles")` to avoid docking the same molecule twice.

---

### Notes / Assumptions
- The receptor used for docking is the **PDBQT** version (`protein_chainA_B_fixed.pdbqt`); PDBFixer produces a fixed **PDB**, but you still must convert it to PDBQT.
- Vina scores are in kcal/mol; they are not directly comparable to `oracle_score` unless you calibrate/validate.
- `ligand_<idx>.pdbqt` filenames depend on the DataFrame index; if you rerun with different filtering/sorting, indices (and filenames) may change.
- The script tracks a “best score” internally, but printing/saving the best entry is currently commented out.
