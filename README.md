
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
- This file is made for Linux environment, you should be able
  to used it on other OS with changing the location of the "download_path"


## USAGE

After cloning this repo onto the Open Enventory server:

1. Change into directory of the new file:
   
   ```bash
   cd update_sql_mol
   ```

2. (Optional): create virtual environment for python to install dependency:
   
   ```bash
   # you can change "update_sql_mol_venv" to other name too
   python3 -m venv update_sql_mol_venv   # Create virtual environment
   source update_sql_mol_venv/bin/activate    # Activate the virtual environment
   ```
  
   if you have **conda**, you can create conda virtual environment with:

   ```bash
   conda create --prefix update_sql_mol_conda-env    # Create virtual environment
   conda activate ./update_sql_mol_conda-env    # Activate the virtual environment
   ```

3. Install python dependencies:
   
   ```bash
   pip install -r requirements.txt
   ```

   if you have used **conda** in step 2 above, you want to do this instead:
   ```bash
   conda install pip -y    # Add pip
   pip install -r requirements.txt    # Install pypi packages
   conda install -c conda-forge molvs -y    # (Optional, for cleaning up mol files) Install molvs using conda, it will automatically install rdkit
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
