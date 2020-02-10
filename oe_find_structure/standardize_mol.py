from rdkit import Chem
from molvs import Standardizer
from pathlib import Path

def standardize_mol(mol_file):

    if Path(mol_file).exists():
        '''Chem.MolFromMolFile() only works with string, not Path object'''
        mol_file = str(mol_file)
        mol = Chem.MolFromMolFile(mol_file)

        s = Standardizer()
        smol = s.standardize(mol)

        with open(mol_file, 'w') as f:
            f.write(Chem.MolToMolBlock(smol))
    
    else:
        # print('file does not exist.')
        raise RuntimeError('File does not exist.')


if __name__ == '__main__':
    pass
    # mol_file = 'data/12259-21-1.mol'
    # standardize_mol(mol_file=mol_file)

