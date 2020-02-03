from rdkit import Chem
from molvs import Standardizer
from pathlib import Path

def standardize_mol(mol_file):
    '''Chem.MolFromMolFile() only works with string, not Path object'''
    mol_file = str(mol_file)
    mol = Chem.MolFromMolFile(mol_file)

    s = Standardizer()
    smol = s.standardize(mol)

    with open(mol_file, 'w') as f:
        f.write(Chem.MolToMolBlock(smol))


if __name__ == '__main__':
    mol_file = 'data/159857-81-5.mol'
    standardize_mol(mol_file=mol_file)
