[![Python 3](https://pyup.io/repos/github/khoivan88/oe_find_structure/python-3-shield.svg)](https://pyup.io/repos/github/khoivan88/oe_find_structure/)
[![Updates](https://pyup.io/repos/github/khoivan88/oe_find_structure/shield.svg)](https://pyup.io/repos/github/khoivan88/oe_find_structure/)

# FIND MISSING STRUCTURE FOR CHEMICALS IN OPEN ENVENTORY

This program is designed specifically for Open Enventory to fix issue with
molecules missing structures (could not be extracted through "Read data from supplier")


## DETAILS

This programs does:

1. Connect into mysql database and find molecule in 'molecule' table
of specific database and find those molecule with missing structure (smiles)
2. folder "missing_mol_files" needs to be created inside /var/lib/mysql
with 'mysql' as ownner (chown mysql:mysql)
3. Try to download mol files from various sources into a folder in
    `/var/lib/mysql/missing_mol_files`
4. Update those sql entries with new downloaded mol_files


## REQUIREMENTS

- root access to the server hosting Open Enventory
- Python 3+
- conda (Optional)
  - conda is used to install `rdkit` and `molvs` to clean mol files (e.g. convert explicit hydrogens to implicit hydrogens, etc.)
  - If you don't already have conda, you can install it using the following <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank">link</a>. If you are not sure what to install, I suggest you install `miniconda3` and **NOT** `anaconda3` for much smaller package footprint.
- This file is made for **Linux** environment, you should be able
  to used it on other OS by changing the location of the ["download_path"](oe_find_structure/find_structure.py#L85)


## USAGE

1. Clone this repository:
   
   ```bash
   git clone https://github.com/khoivan88/oe_find_structure.git    #if you have git
   # if you don't have git, you can download the zip file then unzip
   ```

2. Change into the directory of the program:
   
   ```bash
   cd oe_find_structure
   ```

#### Without conda installed:
   Skip ahead to [this](#with-conda-installed) if you have conda installed.

3. (Optional): create virtual environment for python to install dependency:
   
   ```bash
   # you can change "update_sql_mol_venv" to other name too
   python3 -m venv oe_find_structure_venv   # Create virtual environment
   source oe_find_structure_venv/bin/activate    # Activate the virtual environment on Linux
   # oe_find_structure_venv/Scripts/activate    # Activate the virtual environment on Windows
   ```

4. Install python dependencies:
   
   ```bash
   pip install -r requirements.txt   # Install all dependencies (without rdkit and molvs)
   ```

#### With conda installed:

   Instead of **step 2 AND step 3** above, if you have conda installed, you can do this instead:

   ```bash
   conda env create --prefix oe_find_structure_conda-env --file ./environment.yml    # Create virtual  environment with conda and install all dependancies
   conda activate ./oe_find_structure_conda-env    # Activate the virtual environment
   ```

5. Run the program:
   
   ```bash
   python oe_find_structure/find_structure.py
   ```

   - Answer questions for:
     - mySQL root password (typing password will not be shown on screen)
     - the name of the database you want to update (twice to confirm)
     - url path for your Open Enventory server (including 'http/https' and no trailing '/')

   You can enable debug mode (more error printing during structure search) by adding '`-d`' :
   ```bash
   python oe_find_structure/find_structure.py -d    # Enable debug mode
   ```
