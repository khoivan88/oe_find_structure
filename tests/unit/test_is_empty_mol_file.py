# Add path for the package directory. Otherwise, import will complain.
import sys, os
sys.path.append(os.path.realpath('oe_find_structure'))

from pathlib import Path
import shutil

import pytest
from oe_find_structure.find_structure import is_empty_mol_file


@pytest.fixture
def create_empty_mol_file1():
    '''This file has marker of mol file but contains an empty structure
    Specifically, the content of this file is:
    "
      SciTegic12151716442D

      0  0  0  0  0  0            999 V2000
      M  END
    "
    '''
    Path("tests/tmp").mkdir(parents=True, exist_ok=True)
    with open('tests/tmp/empty_mol_file1.mol', 'w') as file:
        file.write('''        
  SciTegic12151716442D

  0  0  0  0  0  0            999 V2000
  M  END
            ''')


def test_empty_mol_file1(create_empty_mol_file1):
    file = Path('tests/tmp/empty_mol_file1.mol')
    is_empty = is_empty_mol_file(file)
    # Remove file:
    Path.unlink(file)

    assert is_empty == True


@pytest.fixture
def create_empty_mol_file2():
    '''This is a file with just 2 empty lines ('\\n\\n')'''
    Path("tests/tmp").mkdir(parents=True, exist_ok=True)
    with open('tests/tmp/empty_mol_file2.mol', 'w') as file:
        file.write('''        
            ''')


def test_empty_mol_file2(create_empty_mol_file2):
    file = Path('tests/tmp/empty_mol_file2.mol')
    is_empty = is_empty_mol_file(file)
    # Remove file:
    Path.unlink(file)

    assert is_empty == True


@pytest.fixture
def create_normal_mol_file():
    '''This is a normal mol file (CAS#: 159857-81-5).
    This compound has explicit hydrogens drawing'''
    Path("tests/tmp").mkdir(parents=True, exist_ok=True)
    # # print(os.getcwd())
    shutil.copyfile('tests/fixture/159857-81-5-original.mol', 'tests/tmp/159857-81-5.mol')


def test_normal_mol_file(create_normal_mol_file):
    file = Path('tests/tmp/159857-81-5.mol')
    is_empty = is_empty_mol_file(file)
    # Remove file:
    Path.unlink(file)
    assert is_empty == False
