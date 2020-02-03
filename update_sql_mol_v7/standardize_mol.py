from rdkit import Chem
from molvs import Standardizer
from pathlib import Path

def standardize_mol(mol_file):
    
    mol_file = str(mol_file)
    mol = Chem.MolFromMolFile(mol_file)

    s = Standardizer()
    smol = s.standardize(mol)

    with open(mol_file, 'w') as f:
        f.write(Chem.MolToMolBlock(smol))

if __name__ == '__main__':
    # mol_file = 'data/159857-81-5.mol'
    # mol_file = Path('/var/lib/mysql/missing_mol_files' + '/' + '159857-81-5.mol')
    
    mol_file = Path('/home/khoi/programs_for_OE/update_sql_mol/update_sql_mol_v7') / '159857-81-5.mol'
    # print(type(mol_file))

    # import os
    # mol_file = os.path.join('/home/khoi/programs_for_OE/update_sql_mol/update_sql_mol_v7', '159857-81-5.mol')
    # print(type(mol_file))

    # mol_file = '/var/lib/mysql/missing_mol_files/159857-81-5.mol'
    # standardize_mol(mol_file=str(mol_file))
    standardize_mol(mol_file=mol_file)

