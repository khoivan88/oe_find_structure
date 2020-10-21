# Add path for the package directory. Otherwise, import will complain.
import sys, os
sys.path.append(os.path.realpath('oe_find_structure'))

from pathlib import Path

import pytest
from oe_find_structure.find_structure import is_binary_string

@pytest.fixture
def create_normal_string(tmpdir):
    '''This is a normal text string. Content:
    "A normal string.\nAnother very normal string."'''
#     with open('tests/tmp/normal_string.txt', 'w') as file:
    with open(Path(tmpdir) / 'normal_string.txt', 'w') as file:
        file.write('A normal string.\nAnother very normal string.')

def test_normal_string(create_normal_string, tmpdir):
    file = Path(tmpdir) / 'normal_string.txt'
    is_binary = is_binary_string(open(file, 'rb').read(1024))
    # Remove file:
    Path.unlink(file)

    assert is_binary == False



@pytest.fixture
def create_binary_string(tmpdir):
    '''This is a binary string'''
    with open(Path(tmpdir) / 'binary_string', 'wb') as file:
        newFileByteArray = bytearray([123, 3, 255, 0, 100])
        file.write(newFileByteArray)

def test_binary_string(create_binary_string, tmpdir):
    file = Path(tmpdir) / 'binary_string')
    is_binary = is_binary_string(open(file, 'rb').read(1024))
    # Remove file:
    Path.unlink(file)
    assert is_binary == True

