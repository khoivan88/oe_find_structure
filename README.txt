This program is designed specifically for Open Enventory to fix issue with
molecule missing structures (could not be extracted through "Read data from supplier")

This programs does:
    1. Connect into mysql database and find molecule in 'molecule' table
of specific database and find those molecule with missing structure (smiles)
    2. folder "missing_mol_files" needs to be created inside /var/lib/mysql
with 'mysql' as ownner (chown mysql:mysql)
    3. Try to download mol files from chemicalbook.com into a folder in
     /var/lib/mysql/missing_mol_files
    4. Update those sql entries with new downloaded mol mol_files

Note: after running this program, YOU STILL NEED TO LOGIN INTO OPEN ENVENTORY
AS ROOT, GO TO "SETTINGS/BATCH PROCESSING". CHOOSE DATABASE THAT YOU WANT TO Update
AND THEN CHOOSE
    "MOLECULE", "EMPIRICAL FORMULA", "MW", "DEG. OF UNSAT." , "STRUCTURE", AND
    "SMILES"
AND THEN SUBMIT TO UPDATE THE SQL QUERY. ONLY AFTER THIS, THE STRUCTURE WILL SHOW UP

Version 4: changing note:
    - Using simple download bar with bar=wget.bar_thermometer in wget.download
    - Added handling database connection error (wrong password, wrong database)
    - Added confirmation for database
