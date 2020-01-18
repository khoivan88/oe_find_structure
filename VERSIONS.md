
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
