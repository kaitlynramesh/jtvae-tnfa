## README — TNFA (2AZ5) Docking Pipelines: (1) Top-5 Known Binders by *K*ᵢ and (2) Novel Candidates by `oracle_score`

### Overview
These notebooks/scripts implement two related AutoDock Vina docking workflows against TNFA (PDB **2AZ5**, **chains A and B**):

1. **Known binders workflow (Top 5 by experimental *K*ᵢ)**  
   - Reads a TSV of known ligands with *K*ᵢ (nM)  
   - Selects the **5 lowest *K*ᵢ** ligands  
   - Prepares each ligand (SMILES → 3D → PDBQT) and docks with Vina  
   - Separately prepares and docks a provided **native ligand** SMILES

2. **Novel candidates workflow (Top/Bottom by `oracle_score`)**  
   - (Optionally) prepares receptor PDB using **PDBFixer** (adds missing atoms/residues + hydrogens)  
   - Loads a CSV of candidate ligands with `oracle_score`  
   - Selects **top 10** + **bottom 10** by score (deduplicated by SMILES)  
   - Prepares and docks each ligand with Vina

Both workflows use the same docking box and Vina settings.

---

## Inputs

### Common receptor/typing inputs
Under:
- `input_dir = .../final_docking_files/inputs`

Required:
- `ad4_types.json`  
  Atom typing parameters used by Meeko to write ligand PDBQT.
- `protein_chainA_B_fixed.pdbqt`  
  Receptor PDBQT used by AutoDock Vina.

Optional / used for receptor preparation:
- `2az5_chainA_B.pdb`  
  Receptor PDB (chains A/B only) used as the starting point for PDBFixer in the novel-candidates workflow.

### Known binders workflow inputs
- `known_ligands/TNFA_CBB_FINAL_PROJECT.tsv` (tab-separated)  
  Must contain columns:
  - `Ki (nM)`
  - `Ligand SMILES`

Native ligand is provided directly as a SMILES string in the code.

### Novel candidates workflow inputs
- `top_10_novel_candidates.csv`  
  Must contain columns:
  - `smiles`
  - `oracle_score`

---

## Outputs

### Known binders workflow outputs
Under:
- `output_dir = .../final_docking_files/output`

Generated:
- `top5_known_database/ligand_top5_<i>.pdbqt` (prepared ligands)
- `top5_known_database/docked_top5_<i>.pdbqt` (docked poses; top 5 written)
- `top5_known_docking_scores.csv` (rank, smiles, vina_score, time_sec)
- `native/native_ligand.pdbqt`
- `native/native_docked_ligand.pdbqt`
- `native/native_ligand_score.txt`

### Novel candidates workflow outputs
Under:
- `output_dir = .../test_ligands_docking_output`

Generated:
- `ligand_<idx>.pdbqt`
- `docked_<idx>.pdbqt` (top 5 poses written)
- `docking_scores.csv` with:
  - `ligand_id`
  - `ligand_file`
  - `smiles`
  - `vina_score`
  - `time_sec`

---

## Critical step: receptor PDB → PDBQT conversion (command line via Open Babel)

AutoDock Vina requires the receptor in **PDBQT** format. In these workflows, ligand PDBQT files are produced automatically (RDKit + Meeko), but the **receptor PDBQT must be created separately**.

After generating or updating the fixed receptor PDB (e.g., from PDBFixer), run **Open Babel** on the command line:

```bash
obabel protein_fixed.pdb -O protein_fixed.pdbqt -xr
```

Key point:
- The `-xr` flag adds polar hydrogens and merges non-polar hydrogens, which is expected for Vina-style docking.

In your paths, this corresponds to converting:
- `protein_chainA_B_fixed.pdb` → `protein_chainA_B_fixed.pdbqt`

Vina then docks using:
- `protein_pdbqt = .../protein_chainA_B_fixed.pdbqt`

---

## Receptor preparation (used in novel-candidates workflow)
Function: `prepare_receptor(pdb_in, pdb_out)`
- Uses **PDBFixer** to:
  - find missing residues
  - find missing atoms
  - add missing atoms
  - add hydrogens at **pH 7.4**
- Writes `protein_chainA_B_fixed.pdb`

You must still convert that PDB to PDBQT using the **obabel** command above before docking.

---

## Ligand preparation (both workflows)
Function: `smiles_to_pdbqt(smiles, output_pdbqt, atom_params)`
- RDKit:
  - Parse SMILES
  - Add hydrogens
  - 3D embed (ETKDG)
  - UFF optimize geometry
- Meeko:
  - Assign AD4 atom types/rotatable bonds using `ad4_types.json`
  - Write ligand PDBQT

If RDKit cannot parse a SMILES string, the ligand is skipped.

---

## Docking setup (both workflows)
- Vina scoring: `sf_name="vina"`
- Grid maps computed once per run:
  - `center = [-19.28, 74.56, 33.93]`
  - `box_size = [20, 20, 20]`
- Docking parameters:
  - `exhaustiveness = 16`
  - `n_poses = 10` generated
  - top `n_poses = 5` written to disk
- Reported score:
  - best-pose energy from `v.energies(n_poses=1)[0][0]` (kcal/mol; more negative is better, Phillips, R., Kondev, J., Theriot, J., & Orme, N. (2013). Entropy Rules! *Physical Biology of the Cell* (241-244). Garland Science. https://books.google.com/books?id=JnyPZwEACAAJ)


---

## Selection logic

### Known binders workflow (Top 5 by *K*ᵢ)
1. Read TSV (`sep="\t"`)
2. Sort by `Ki (nM)` ascending
3. Take top 5 `Ligand SMILES`
4. Dock each ligand and save scores
5. Separately prepare + dock the provided native ligand SMILES

### Novel candidates workflow (Top/Bottom by `oracle_score`)
1. Read CSV
2. Coerce `oracle_score` to numeric; drop missing
3. Sort descending by `oracle_score`
4. Take top 10 and bottom 10, concatenate
5. Drop duplicate SMILES
6. Dock each ligand and save results to `docking_scores.csv`

---

## Notes / assumptions
- Vina scores (kcal/mol) are not directly comparable to *K*ᵢ or your `oracle_score` without validation/calibration.
- Receptor correctness (protonation, missing loops, chain selection) strongly affects docking; PDBFixer helps but does not guarantee a biologically perfect model.
- Filenames like `ligand_<idx>.pdbqt` depend on the DataFrame index; reruns can change indices unless stabilized.
