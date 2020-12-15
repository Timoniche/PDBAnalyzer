from Bio.PDB.PDBParser import PDBParser

pdb_reader = PDBParser(PERMISSIVE=1)
structure_id = "3DMax_chr1.pdb"
filename = "pdb_files/3DMax_chr1.pdb"
structure = pdb_reader.get_structure(structure_id, filename)

model = structure[0]

for chain in model:
    for residue in chain:
        for atom in residue:
            if atom.get_name() == 'CA':
                print(atom.get_id(), residue.get_resname(), residue.get_id()[1], chain.get_id(), atom.get_coord())
