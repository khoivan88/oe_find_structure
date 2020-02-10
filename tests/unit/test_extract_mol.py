# Add path for the package directory. Otherwise, import will complain.
import sys, os
sys.path.append(os.path.realpath('oe_find_structure'))

import pytest
from oe_find_structure.find_structure import extract_mol_from_pubchem, \
                                            extract_mol_from_cactus, \
                                            extract_mol_from_chemicalbook, \
                                            extract_mol


@pytest.mark.parametrize(
    "cas, return_value", [
        ('159857-81-5', 0),
        ('12259-21-1', 0),
    ]
)
def test_extract_mol_from_pubchem_without_existing_mol_files(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_pubchem() WITHOUT existing mol files'''

    '''Changing the value of 'download_path' variable to tmpdir'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)

    result = extract_mol_from_pubchem(cas)  
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('159857-81-5', -1),
        ('12259-21-1', -1),
    ]
)
def test_extract_mol_from_pubchem_with_existing_mol_files(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_pubchem() WITH existing mol files'''
    
    '''Changing the value of 'download_path' variable'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)
    
    # '''Create these files first to generate condition of existing files'''
    extract_mol_from_pubchem(cas)

    result = extract_mol_from_pubchem(cas)
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('00000-00-0', '00000-00-0'),
        ('14636-75-7', '14636-75-7'),
    ]
)
def test_extract_mol_from_pubchem_not_found(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_pubchem() for chemicals NOT exist in Pubchem'''

    '''Changing the value of 'download_path' variable'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)

    result = extract_mol_from_pubchem(cas)
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('53597-25-4', '53597-25-4'),
        ('9035-81-8', '9035-81-8'),
        ('70750-57-1', '70750-57-1'),
    ]
)
def test_extract_mol_from_pubchem_found_but_no_structure(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_pubchem() for chemicals found in Pubchem but without structures'''

    '''Changing the value of 'download_path' variable'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)

    result = extract_mol_from_pubchem(cas)
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('12259-21-1', 0),
        ('10170-69-1', 0)
    ]
)
def test_extract_mol_from_cactus_without_existing_mol_files(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_cactus() WITHOUT existing mol files'''

    '''Changing the value of 'download_path' variable to tmpdir'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)

    result = extract_mol_from_cactus(cas)  
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('12259-21-1', -1),
        ('10170-69-1', -1)
    ]
)
def test_extract_mol_from_cactus_with_existing_mol_files(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_cactus() WITH existing mol files'''
    
    '''Changing the value of 'download_path' variable'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)
    
    # '''Create these files first to generate condition of existing files'''
    extract_mol_from_cactus(cas)

    result = extract_mol_from_cactus(cas)
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('00000-00-0', '00000-00-0'),
        ('12073-36-8', '12073-36-8'),
    ]
)
def test_extract_mol_from_cactus_not_found(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_cactus() for chemicals NOT exist in Pubchem'''

    '''Changing the value of 'download_path' variable'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)

    result = extract_mol_from_cactus(cas)
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('159857-81-5', 0),
        ('12073-36-8', 0),
        ('10170-69-1', 0)
    ]
)
def test_extract_mol_from_chemicalbook_without_existing_mol_files(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_chemicalbook() WITHOUT existing mol files'''

    '''Changing the value of 'download_path' variable to tmpdir'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)

    result = extract_mol_from_chemicalbook(cas)  
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('159857-81-5', -1),
        ('12073-36-8', -1),
        ('10170-69-1', -1)
    ]
)
def test_extract_mol_from_chemicalbook_with_existing_mol_files(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_chemicalbook() WITH existing mol files'''
    
    '''Changing the value of 'download_path' variable'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)
    
    # '''Create these files first to generate condition of existing files'''
    extract_mol_from_chemicalbook(cas)

    result = extract_mol_from_chemicalbook(cas)
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('00000-00-0', '00000-00-0'),
        ('14636-75-7', '14636-75-7'),
    ]
)
def test_extract_mol_from_chemicalbook_not_found(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_chemicalbook() for chemicals NOT exist in Pubchem'''

    '''Changing the value of 'download_path' variable'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)

    result = extract_mol_from_chemicalbook(cas)
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('12259-21-1', '12259-21-1'),
        ('53597-25-4', '53597-25-4'),
        ('70750-57-1', '70750-57-1'),
    ]
)
def test_extract_mol_from_chemicalbook_found_but_no_structure(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_chemicalbook() for chemicals found in chemicalbook but without structures'''

    '''Changing the value of 'download_path' variable'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir)

    result = extract_mol_from_chemicalbook(cas)
    assert result == return_value


@pytest.mark.parametrize(
    "cas, return_value", [
        ('72773-04-7', 0),    # found in chemicalbook
        ('159857-81-5', 0),    # found in chemicalbook
        ('12259-21-1', 0),    # found in cactus
        ('9005-84-9', 0),     # found in cactus
        # ('870987-63-6', 0),    # found in pubchem
        ('53597-25-4', '53597-25-4'),    # structure not found in any database
        ('70750-57-1', '70750-57-1'),    # structure not found in any database
        ('00000-00-0', '00000-00-0')
    ]
)
def test_extract_mol_combined(tmpdir, monkeypatch, cas, return_value):
    '''Test extract_mol_from_chemicalbook() for chemicals found in chemicalbook but without structures'''

    '''Changing the value of 'download_path' variable'''
    # print('tmpdir is: {}'.format(tmpdir))
    monkeypatch.setattr("oe_find_structure.find_structure.download_path", tmpdir + '/')

    result = extract_mol(cas)
    assert result == return_value
