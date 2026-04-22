# Physical metrics for validating TNFa binders
Hello Reader! Welcome to my section of the code. Note that most of this relies heavily on the Vina evaluations done before this, so be sure to check that out first. In terms of finding del-rSASA, the provided .ipynb file should provide you with the essential components in a standard conda baseenv. The contacts were done in a bit more tricky of a manner and are described below.

To obtains contacts, open solo_protein and any of the ligands in PyMol. Split the ligand states in the GUI, followed by deleting all but the state of interest (typically the first). Do the following set of commands:

select chain A or docking_X_000Y // X and Y being the ligand of interest, and state of interest respectively

In the selection window and 'action' GUI, 'find' contacts within 4.0 A. Rename the contact map produced 'A_ligand_contacts.'

Redo the above with chain B.

Finally run the following set of commands:

x = cmd.get_session('A_ligand_contacts', 1)['names'][0]
print(x[-2][-2][0][8])

At this point pull the subsequent list as given, and find the length. (Suffices to put into len() in provided .ipynb) This is the number of contacts with chain A.

x = cmd.get_session('B_ligand_contacts', 1)['names'][0]
print(x[-2][-2][0][8])

As before, copy the list, find its length and record as the number of contacts with chain B. Add the two results and obtain the number of contacts.
