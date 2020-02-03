
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
     /var/lib/mysql/missing_mol_files
4. Update those sql entries with new downloaded mol_files


## REQUIREMENTS

- root access to the server hosting Open Enventory
- Python 3+
- conda (Optional)
  - conda is used to install `rdkit` and `molvs` to clean mol files (e.g. convert explicit hydrogens to implicit hydrogens, etc.)
  - If you don't already have conda, you can install it using the following <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank">link</a>. If you are not sure what to install, I suggest you install `miniconda3` and **NOT** `anaconda3` for much smaller package footprint.
- This file is made for **Linux** environment, you should be able
  to used it on other OS by changing the location of the ["download_path"](update_sql_mol_v7/update_sql_mol.py#L88)



## USAGE

After cloning this repo onto the Open Enventory server:

1. Change into the directory of the program:
   
   ```bash
   cd update_sql_mol
   ```

#### Without conda installed:
   Skip ahead to [this](#with-conda-installed) if you have conda installed.

2. (Optional): create virtual environment for python to install dependency:
   
   ```bash
   # you can change "update_sql_mol_venv" to other name too
   python3 -m venv update_sql_mol_venv   # Create virtual environment
   source update_sql_mol_venv/bin/activate    # Activate the virtual environment
   ```

3. Install python dependencies:
   
   ```bash
   pip install -r requirements.txt   # Install all dependencies (without rdkit and molvs)
   ```

#### With conda installed:

   Instead of **step 2 AND step 3** above, if you have conda installed, you can do this instead:

   ```bash
   conda env create --prefix update_sql_mol_conda-env --file ./environment.yml    # Create virtual  environment with conda and install all dependancies
   conda activate ./update_sql_mol_conda-env    # Activate the virtual environment
   ```

4. Run the program:
   
   ```bash
   python update_sql_mol_v7/update_sql_mol.py    # Replace "update_sql_mol_v6" with latest version
   ```

   - Answer questions for:
     - confirming running under root
     - mySQL root password (typing password will not be shown on screen)
     - the name of the database you want to update (twice to confirm)
     - url path for your Open Enventory server (including 'http/https' and no trailing '/')
