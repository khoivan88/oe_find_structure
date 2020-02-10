
## Version 4:

- Using simple download bar with bar=wget.bar_thermometer in wget.download
- Added handling database connection error (wrong password, wrong database)
- Added confirmation for database


## Version 5

- Change to https (from http) for ChemicalBook site


## Version 6 (2020-01-18):

- Remove this extra step of going into Open Enventory manual by the user
using a web browser and doing Batch Processing.

  Why: for versions 5 and before:
>       after running this program,
>       YOU STILL NEED TO LOGIN INTO OPEN ENVENTORY ON A WEB BROWSER AS ROOT,
>       GO TO "SETTINGS/BATCH PROCESSING". CHOOSE DATABASE THAT YOU WANT TO UPDATE
>       AND THEN CHOOSE
>           "MOLECULE", "EMPIRICAL FORMULA", "MW", "DEG. OF UNSAT." ,
>           "STRUCTURE", "SMILES" AND "FINGERPRINT"
>       AND THEN SUBMIT TO UPDATE THE SQL QUERY.
>       ONLY AFTER THIS, THE STRUCTURE WILL SHOW UP


## Version 7 (2020-02-02):

- Add some check to make sure the mol file is valid (not a binary string or empty mol files)
- Add searching for SD File (similar to mol file) from Cactus: https://cactus.nci.nih.gov/chemical/structure
- Add searching for SD File (similar to mol file) from pubchem: https://pubchem.ncbi.nlm.nih.gov/
- Add function to clean mol file (e.g: remove explicit hydrogens, etc.). User needs to have rdkit and molvs libraries installed first. Suggest install these libraries using conda. If these libraries are not found, this program will skip this function and use the mol files as is.
- Run looking for mol files function twice. The first time with Pool(20) to take advantage of multiple threads for sites without request limit. The second time, set to 'Pool()' to take advantage of Pubchem, a large collection but has limited (no more than 5 requests per second) request rate.
- Reduce error output and add debug mode (print more error)

