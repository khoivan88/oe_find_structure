import pytest
import shutil
# import os
# from setuptools import setup, find_packages
from oe_find_structure.standardize_mol import standardize_mol
from pathlib import Path

@pytest.fixture
def create_bad_file(tmpdir):
    '''This is a SDF file with a lot more info than mol file 
    and it has explicit hydrogens drawing'''
    new_file = Path(tmpdir) / '159857-81-5.mol'

    # make a duplicate to preserve file
    shutil.copyfile('tests/fixture/159857-81-5-original.mol', new_file)
    return new_file


@pytest.fixture
def create_good_file(tmpdir):
    '''This is a mol file with implicit hydrogens drawing'''
    new_file = Path(tmpdir) / '159857-81-5-good.mol'

    # make a duplicate to preserve file
    shutil.copyfile('tests/fixture/159857-81-5-implicitH.mol', new_file)
    return new_file


def test_file_not_found(tmpdir):
    file = Path(tmpdir) / 'non-existing.mol'
    # print(file)
    with pytest.raises(RuntimeError) as error:
        standardize_mol(file)
    assert 'File does not exist.' in str(error.value)


def test_remove_explicitH(create_bad_file):
    correct_mol_content = open('tests/fixture/159857-81-5-implicitH.mol', 'r').read()
    
    file = Path(create_bad_file)
    # print(file)
    assert open(file, 'r').read() != correct_mol_content

    standardize_mol(file)
    mol_content_fixed = open(file, 'r').read()

    assert mol_content_fixed == correct_mol_content


def test_no_change(create_good_file):
    correct_mol_content = open('tests/fixture/159857-81-5-implicitH.mol', 'r').read()

    file = Path(create_good_file)
    assert open(file, 'r').read() == correct_mol_content

    standardize_mol(file)
    mol_content_fixed = open(file, 'r').read()

    assert mol_content_fixed == correct_mol_content
