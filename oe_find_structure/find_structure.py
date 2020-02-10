#!/usr/bin/python

"""
This program is designed specifically for Open Enventory to fix issue with
molecule missing structures (could not be extracted through "Read data from supplier")

This programs does:
    1. Connect into mysql database and find molecule in 'molecule' table
of specific database and find those molecule with missing structure (smiles)
    2. folder "missing_mol_files" needs to be created inside /var/lib/mysql
with 'mysql' as ownner (chown mysql:mysql)
    3. Try to download mol files from chemicalbook.com into a folder in
     /var/lib/mysql/missing_mol_files
    4. Update those sql entries with new downloaded mol_files

Update for version 6: this extra step of going into Open Enventory using a web
browser and doing Batch Processing is not neccessary anymore
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

Version 5: changing note:
    - Change to https (from http) for ChemicalBook site

Version 6: changing note:
    - Remove this extra step of going into Open Enventory manual by the user 
    using a web browser and doing Batch Processing.
    Why: for Version 5 and before: after running this program, 
    >     YOU STILL NEED TO LOGIN INTO OPEN ENVENTORY ON A WEB BROWSER AS ROOT, 
    >     GO TO "SETTINGS/BATCH PROCESSING". CHOOSE DATABASE THAT YOU WANT TO UPDATE
    >     AND THEN CHOOSE
    >         "MOLECULE", "EMPIRICAL FORMULA", "MW", "DEG. OF UNSAT." , 
    >         "STRUCTURE", "SMILES" AND "FINGERPRINT"
    >     AND THEN SUBMIT TO UPDATE THE SQL QUERY. 
    >     ONLY AFTER THIS, THE STRUCTURE WILL SHOW UP

Version 7:
    - Add some check to make sure the mol file is valid 
    (not a binary string or empty mol files)
    - Add searching for SD File (similar to mol file) from Cactus:
    https://cactus.nci.nih.gov/chemical/structure
    - Add searching for SD File (similar to mol file) from pubchem:
    https://pubchem.ncbi.nlm.nih.gov/
    - Add function to clean mol file (e.g: remove explicit hydrogens, etc.).
    User needs to have rdkit and molvs libraries installed first. 
    Suggest install these libraries using conda. If these libraries are not found, 
    this program will skip this function and use the mol files as is.
    - Run looking for mol files function twice. The first time with Pool(20)
    to take advantage of multiple threads for sites without request limit.
    The second time, set to 'Pool()' to take advantage of Pubchem, a large 
    collection but has limited (no more than 5 requests per second) request rate.
"""

import getpass
import importlib
import itertools
import os
import re
import sys
import threading
import time
from multiprocessing import Pool
from pathlib import Path

import mysql.connector as mariadb
import pubchempy as pcp  # https://pubchempy.readthedocs.io/en/latest/guide/gettingstarted.html
import requests

# Check if rkdit and molvs libraries are available
rdkit_spec = importlib.util.find_spec("rdkit")
molvs_spec = importlib.util.find_spec("molvs")
lib_found = (rdkit_spec is not None) and (molvs_spec is not None)
if lib_found:
    from standardize_mol import standardize_mol

missing_mol_file = set()
download_path = '/var/lib/mysql/missing_mol_files'


def main():
    global download_path

    # Require user running this python as root
    if_root = input('Are you login as root user? (y/n): ')
    if (if_root not in ['y', 'yes']):
        print('You need to convert to root user before running this program')
        exit(1)

    # Get user input for root password and the database needs to be updated
    # to hide password input: https://stackoverflow.com/questions/9202224/getting-command-line-password-input-in-python
    # password = input('Please type in the password for "root" user: ')
    password = getpass.getpass('Please type in the password for MySQL "root" user: ')
    database = input('Please type in the name of the database needs updating: ')
    # Ask user to retype the database name and if it does NOT match, exit the programs
    database2 = input('Please re-type the name of the database to confirm: ')
    if (database != database2):
        print('Database names do NOT match!')
        exit(2)
    # Ask user for OE URL path
    oe_url_path_input = input(
        '''Please type the URL address of your database (including http/https)
        (leave blank if not sure, the default will be: 'localhost'): ''')

    """
    Info for mysql connection and query can be found here:
     https://mariadb.com/resources/blog/how-connect-python-programs-mariadb
     https://dev.mysql.com/doc/connector-python/en/connector-python-tutorial-cursorbuffered.html

    Handling error in password or database not exists:
    https://dev.mysql.com/doc/connector-python/en/connector-python-example-connecting.html
    https://dev.mysql.com/doc/connector-python/en/connector-python-api-errorcode.html
    """

    # Open a connection to mysql
    try:
        mariadb_connection = mariadb.connect(user='root', password=password, database=database)
        # Create a cursor in the sql table using the open connection
        cursor_select = mariadb_connection.cursor(buffered=True)

        '''-----Step1: run SELECT query to find CAS# for those with missing structures-----'''
        print('Getting molecule with missing structures. Please wait!')
        # query = ("SELECT distinct cas_nr FROM molecule WHERE smiles='' and cas_nr!=''")
        # query = ("SELECT distinct cas_nr FROM molecule WHERE cas_nr!='' and molfile_blob like '%open enventory%'")
        query = ("SELECT distinct cas_nr FROM molecule WHERE cas_nr!='' and smiles=''")
        try:
            cursor_select.execute(query)
        except mariadb.Error as error:
            print('Error: {}'.format(error))

        # Create a new empty list and added the cas# into this new list:
        to_be_downloaded = []
        for (cas_nr, ) in cursor_select:
            to_be_downloaded.append(cas_nr)

        '''-----Step 2: downloading mol file-----'''
        # Check if download path with the missing_mol_file directory exists. If not, create it
        # https://stackoverflow.com/questions/12517451/automatically-creating-directories-with-file-output
        # https://docs.python.org/3/library/os.html#os.makedirs
        os.makedirs(download_path, exist_ok=True)

        print('Downloading missing mol files. Please wait!')
        try:
            '''
            For the first run, you can set Pool(10) or Pool(20) to speed up process,
            however, after the first 1 or two run, set Pool() to take advantage of Pubchem
            Reason: Pubchem blocks access if there are more than 5 request per seconds
            '''
            with Pool(20) as p:
                still_missing_first_round = p.map(extract_mol, to_be_downloaded)
                # 'still_missing_list' is a list of the return value from extract_mol().
                # This function return CAS (string) of chemicals whose mol file cannot be found
            
            still_missing_first_round = [x for x in still_missing_first_round if x]
            # print('Still missing chemicals after first round: {}'.format(still_missing_first_round))

            '''Run the extract_mol() the second time with Pool() to be 
            able to take advantage of Pubchem service. Pubchem has a large collection.
            However, Pubchem API has a limit of not more than 5 request per seconds'''
            time.sleep(3)
            with Pool() as p:
                still_missing_second_round = p.map(extract_mol, still_missing_first_round)
                # 'still_missing_list' is a list of the return value from extract_mol().
                # This function return CAS (string) of chemicals whose mol file cannot be found

            still_missing_second_round = [x for x in still_missing_second_round if x]
            # print('Still missing chemicals after 2nd round: {}'.format(still_missing_second_round))

        except Exception as error:
            print(error)


        finally:
            '''-----Step 3: run UPDATE query to upload-----'''
            '''Upload mol files into MySQL'''
            print('Updating SQL table!')
            count_file_updated = 0
            for cas in to_be_downloaded:
                try:
                    # run update_sql() and also increment the count for successful update
                    # update_sql() return 1 if successs, otherwise return 0
                    count_file_updated += update_sql(mariadb_connection, cas)

                except mariadb.Error as error:
                    print('Error: {}'.format(error))

            mariadb_connection.close()
            print('\nStill missing mol files:\n{}'.format(missing_mol_file))
            print('\nSummary: ')
            print('\t{} mol files are still missing.'.format(len(missing_mol_file)))
            # path, dirs, files = next(os.walk(download_path))
            # file_count = len(files)
            print('\t{} mol files updated! '.format(count_file_updated))


            '''-----Step 4: Use OE Batch Processing to generate image for structure-----'''
            # Login into OE and use Batch Processing to generate structure images
            oe_url_path = 'http://localhost' if oe_url_path_input == '' else oe_url_path_input
            # Only run the batch processing if there is any updated mol files:
            if count_file_updated > 0: 
                oe_batch_process(
                    database=database,
                    user='root',
                    password=password,
                    oe_url_path=oe_url_path)

            print()
            if not lib_found:
                unfound_lib_string = ("rdkit and molvs libraries not found. " 
                    "mol files are used as is!")
                print(unfound_lib_string)

    except mariadb.Error as err:
        if err.errno == mariadb.errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == mariadb.errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    else:
        mariadb_connection.close()


def update_sql(mariadb_connection, cas_nr):
    global download_path
    cursor_update = mariadb_connection.cursor(buffered=True)
    file_path = download_path + '/{}.mol'.format(cas_nr)
    mol_file = Path(file_path)
    # print(mol_file)

    # result = extract_mol(cas_nr)
    # if molfile exists or downloaded (extract_mol return -1 or 0)
    if mol_file.exists():
        print('CAS# {:<15}: '.format(cas_nr), end='')
        # cursor_update.execute(insert_mol_file, (file_path, cas_nr))
        cursor_update.execute("UPDATE molecule SET molfile_blob=LOAD_FILE('{}') WHERE cas_nr='{}'".format(mol_file, cas_nr))
        mariadb_connection.commit()
        print('mol file uploaded successfully!')
        return 1
    
    # extract_mol return the cas# of those that it could not find mol file
    else:
        missing_mol_file.add(cas_nr)
        return 0


def extract_mol(cas_nr):
    print('Looking for {} ...'.format(cas_nr))

    '''Assume mol file existed until exhaust all searches'''
    mol_file_existed = True

    still_missing_cas = ''

    '''Find mol file from chemicalbook'''
    chemicalbook_result = extract_mol_from_chemicalbook(cas_nr)
    '''cas# (string) is return if mol file cannot be found'''
    if isinstance(chemicalbook_result, str):
        # print('CAS {} not found from chemicalbook.com. Trying cactus now.'.format(cas_nr))
        '''Find mol file from cactus'''
        cactus_result = extract_mol_from_cactus(cas_nr)
        if isinstance(cactus_result, str):
            # print('CAS {} not found from cactus service. '.format(cas_nr))
            '''Find mol file from pubchem'''
            pubchem_result = extract_mol_from_pubchem(cas_nr)
            if isinstance(pubchem_result, str):
                # print('CAS {} not found from PubChem service. '.format(cas_nr))
                '''Exhaust all search so change mol_file_existed to False'''
                mol_file_existed = False
                '''set still_missing_cas'''
                still_missing_cas = pubchem_result


    '''Clean up mol file'''
    global download_path
    global lib_found

    if lib_found and mol_file_existed:
        file_name = cas_nr + '.mol'
        download_file = Path(download_path) / file_name
        try:
            '''Clean up mol file'''
            if download_file.exists() and os.stat(download_file).st_size != 0:
                # print(download_file) 
                # print('Cleaning mol file for {}'.format(download_file))
                standardize_mol(mol_file=download_file)
                return 0

        except Exception as error:
            # print('\nCleaning mol file for {}'.format(download_file))
            print('\tError for {}: {}\n'.format(cas_nr, error))
            # return cas_nr

    '''Return still_missing_cas if exist'''
    # if still_missing_cas is not None and still_missing_cas != '':
    if still_missing_cas:
        return still_missing_cas


def extract_mol_from_chemicalbook(cas_nr):
    global download_path
    '''
    This function is used to extract a single mol file from chemicalbook.com
    See here for more info: http://stackabuse.com/download-files-with-python/
    '''

    headers = {
        'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/53.0.2785.143 Safari/537.36'}

    # get url from ChemicalBook to download mol file
    url = 'https://www.chemicalbook.com/CAS/mol/'
    cas = cas_nr
    file_name = cas + '.mol'
    full_url = url + file_name

    # download_path = '/Users/khoivan/Downloads/mol_files/'
    download_file = Path(download_path) / file_name

    # Check if the file not exists and download
    # check file exists: https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists
    if download_file.exists() and os.stat(download_file).st_size != 0:
        # print('{} already downloaded'.format(file_name))
        return -1
    else:
        # requests.get().history would show code
        # 302 (redirect) if the mol file does not exist.
        # ChemicalBook automatically redirect if could not find mol file
        # redirected = False;
        try:
            result = requests.get(full_url, headers=headers, timeout=10)
            hist_len = len(result.history)
            
            if result.status_code == 200 and hist_len == 0:  # requests.get.history will a list of more than 0 element if redirected
                # print(result.text, file=open('test.html', 'w'))
                download_file.write_text(data=result.text)
                # Check if the mol file is a binary string (some error during downloading) or empty mol file:
                if is_binary_string(open(download_file, 'rb').read(1024)) or is_empty_mol_file(download_file):
                    os.remove(download_file)    # remove the error mol file
                    return cas_nr
                else:
                    return 0

            elif result.status_code == 200 and hist_len > 0:  # requests.get.history will a list of more than 0 element if redirected
                # print(result.text, file=open('test.html', 'w'))
                chemicalbook_search_url = 'https://www.chemicalbook.com/Search_EN.aspx?keyword={}'.format(cas_nr)
                redirect_search = requests.get(chemicalbook_search_url, headers=headers, timeout=10)
                
                if redirect_search.status_code == 200 and len(redirect_search.history) == 0:  # requests.get.history will a list of more than 0 element if redirected
                    mol_file_pattern = re.compile(r'href=\'(.+)\'>Mol file')
                    mol_file_link_suffix = mol_file_pattern.search(redirect_search.text).group(1)
                    # print(mol_file_link_suffix)

                    cas_number_pattern = re.compile(r'(?<![\d\w])(\d{2,7}-\d{2}-\d)')
                    """Explain above regex:
                        "(?<![\d\w])(\d{2,7}-\d{2}-\d)"
                        Negative Lookbehind (?<![\d\w])
                        Assert that the Regex below does not match
                            Match a single character present in the list below [\d\w]
                            \d matches a digit (equal to [0-9])
                            \w matches any word character (equal to [a-zA-Z0-9_])
                        1st Capturing Group (\d{2,7}-\d{2}-\d)
                            \d{2,7} matches a digit (equal to [0-9])
                            {2,7} Quantifier — Matches between 2 and 7 times, as many times as possible, giving back as needed (greedy)
                            - matches the character - literally (case sensitive)
                            \d{2} matches a digit (equal to [0-9])
                            {2} Quantifier — Matches exactly 2 times
                            - matches the character - literally (case sensitive)
                            \d matches a digit (equal to [0-9])
                    """
                    redirect_search_cas = cas_number_pattern.search(mol_file_link_suffix).group(1)
                    # print(redirect_search_cas)

                    # Make sure that the result found by chemicalbook redirect search match the input CAS#
                    if redirect_search_cas == cas_nr:
                        new_mol_file_url = 'https://www.chemicalbook.com/{}'.format(mol_file_link_suffix)
                        new_result = requests.get(new_mol_file_url, headers=headers, timeout=10)
                
                        if new_result.status_code == 200 and len(new_result.history) == 0:  # requests.get.history will a list of more than 0 element if redirected
                            # print(new_result.text)
                            download_file.write_text(data=new_result.text)
                            # Check if the mol file is a binary string (some error during downloading) or empty mol file:
                            if is_binary_string(open(download_file, 'rb').read(1024)) or is_empty_mol_file(download_file):
                                os.remove(download_file)    # remove the error mol file
                                return cas_nr
                            else:
                                return 0
            
            # return the cas # of the chemical whose mol file cannot be found
            return cas

        except Exception as error:
            print(error)
            return cas


def extract_mol_from_cactus(cas_nr):
    global download_path
    '''
    This function is used to extract a single mol file
    from: https://cactus.nci.nih.gov/chemical/structure
    '''

    headers = {
        'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/53.0.2785.143 Safari/537.36'}

    # get url from cactuc to download mol file
    url = 'https://cactus.nci.nih.gov/chemical/structure/{}/file?format=sdf'.format(cas_nr)
    file_name = cas_nr + '.mol'
    download_file = Path(download_path) / file_name
    # print(download_file)

    # Check if the file not exists and download
    # check file exists: https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists
    if download_file.exists() and os.stat(download_file).st_size != 0:
        # print('{} already downloaded'.format(file_name))
        return -1
    else:
        try:
            # print('\tSearching Cactus (NIH) service...')
            result = requests.get(url, headers=headers, timeout=30)
            hist_len = len(result.history)

            # print(result.text)
            if result.status_code == 200 and hist_len == 0:  # requests.get.history will a list of more than 0 element if redirected
                download_file.write_text(data=result.text)
                
                # Check if the mol file is a binary string (some error during downloading) or empty mol file:
                if is_binary_string(open(download_file, 'rb').read(1024)) or is_empty_mol_file(download_file):
                    os.remove(download_file)    # remove the error mol file
                    return cas_nr

                return 0

            # return the cas # of the chemical whose mol file cannot be found
            return cas_nr

        except Exception as error:
            print(error)
            return cas_nr


def extract_mol_from_pubchem(cas_nr):
    global download_path
    headers = {
        'user-agent': 'Mozilla/5.0 (X11; CentOS; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/73.0.3683.75 Safari/537.36'}

    try:
        # print('\tSearching Pubchem...')

        # Using pubchem api for python
        # Getting CID number, the result of this, by default is exact match. The result is returned as a list.
        # cid = pcp.get_cids(cas_nr, 'name', 'substance', list_return='flat')
        cid = pcp.get_cids(cas_nr, 'name')

        file_name = cas_nr + '.mol'
        download_file = Path(download_path) / file_name

        # Check if the file not exists and download
        # check file exists: https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists
        if download_file.exists() and os.stat(download_file).st_size != 0:
            # print('{} already downloaded'.format(file_name))
            return -1
        
        else:

            #  this api return an empty list if it cannot find cas_nr. This is to check if pubchem has this chemical.
            if len(cid) > 0:
                # if Pubchem found the result, get the first result of the list
                cid = cid[0]
                # print('Compound ID (CID) from PubChem is: {} and type is: {}'.format(cid, type(cid)))
            
                # To double check if the CAS number is correct:
                # using pubchem api, get a list of synonym. The result is a list of dict.
                # choose the first result and check first 5 values for 'Synonym' key:
                # synonyms = pcp.get_synonyms(cid)[0]['Synonym'][:7]
                synonyms = pcp.get_synonyms(cid)[0]['Synonym']
                # print('List of synonyms is: {}'.format(synonyms)); exit(0)

                if cas_nr not in synonyms:
                    raise ValueError('\tThis is not an exact match!')

                # get url from Fisher to get url to download sds file
                get_sdf_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/sdf'.format(cid)

                # # Check if the file not exists and download
                # # check file exists: https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists
                # if download_file.exists():
                #     # print('{} already downloaded'.format(file_name))
                #     return -1
                # else:

                # # Another way to get sdf, from pubchempy ---------------------------------------
                #     sdf = pcp.get_sdf(cid)
                #     with open('159857-81-5.mol', 'w') as f:
                #         f.write(sdf)
                # # ----------------------------------------------------------------------------------

                # Get the html request info using CID number from pubchem
                r = requests.get(get_sdf_url, headers=headers, timeout=15)
                # print('url is: {}'.format(get_sdf_url))

                # Check to see if give OK status (200) and not redirect
                if r.status_code == 200 and len(r.history) == 0:
                    download_file.write_text(data=r.text)
                    
                    # Check if the mol file is a binary string (some error during downloading) or empty mol file:
                    if is_binary_string(open(download_file, 'rb').read(1024)) or is_empty_mol_file(download_file):
                        os.remove(download_file)    # remove the error mol file
                        return cas_nr
                    else:
                        return 0
                
            # If not, try to find substances as well
            elif len(cid) == 0:
                '''pcp.get_substances(cas_nr, 'name') returns a list of Substances if found: 
                Ref: https://github.com/mcs07/PubChemPy/blob/e3c4f4a9b6120433e5cc3383464c7a79e9b2b86e/pubchempy.py#L328'''
                substances = pcp.get_substances(cas_nr, 'name')   
                # print(sid); exit(0)

                if len(substances) == 0:
                    # print('nothing here')
                    raise ValueError('Could not find any compounds or substances with this CAS {} on Pubchem.'.format(cas_nr))
                else:
                    for substance in substances:
                        # print('Substance ID (SID) from PubChem is: {} and type is: {}'.format(substance, type(substance)))

                        '''Ref: https://github.com/mcs07/PubChemPy/blob/e3c4f4a9b6120433e5cc3383464c7a79e9b2b86e/pubchempy.py#L735'''
                        # substance_synonyms = substance.to_dict(properties=['synonyms'])['synonyms']
                        '''
                        substance.to_dict(properties=['synonyms']) return example:
                        {'synonyms': ['12259-21-1', 'Iron oxide (Fe2O3), hydrate', 'Ferric oxide hydrate', 
                                        'Ferrox', 'Hydrated ferric oxide', 'Hydrous ferric oxide', 
                                        'Iron oxide (Fe203), hydrate']}
                        '''
                        
                        substance_synonyms = substance.synonyms   # https://github.com/mcs07/PubChemPy/blob/e3c4f4a9b6120433e5cc3383464c7a79e9b2b86e/pubchempy.py#L1095
                        '''
                        substance.synonyms' return example:
                            ['12259-21-1', 'Iron oxide (Fe2O3), hydrate', 'Ferric oxide hydrate', 
                            'Ferrox', 'Hydrated ferric oxide', 'Hydrous ferric oxide', 
                            'Iron oxide (Fe203), hydrate']
                        '''

                        # Check to make sure the substance has the same CAS#
                        if cas_nr in substance_synonyms:
                            sdf = pcp.get_sdf(identifier=substance.sid, namespace='sid', domain='substance')
                            # print(sdf)
                            if sdf:    # pcp.get_sdf return None if not found SDF                               
                                download_file.write_text(data=sdf)
                                
                                # Check if the mol file is a binary string (some error during downloading) or empty mol file:
                                if is_binary_string(open(download_file, 'rb').read(1024)) or is_empty_mol_file(download_file):
                                    os.remove(download_file)    # remove the error mol file
                                else:
                                    return 0

            # If none of the Substances has the same CAS and/or has SDF (mol) file, then return the CAS #
            return cas_nr

    except Exception as error:
        # print('.', end='')
        print('\t{}'.format(error))
        return cas_nr


def oe_batch_process(database: str, password: str, user: str = 'root', oe_url_path: str = 'http://localhost'):
    """This function is used to tell Open Enventory to generate pictures
    of chemical structures.
    After having the mol files in mySQL, the structure of chemicals 
    still need to be generated (into .png file) before they can be seen 
    from OE interface. This function handles that
    Args:
        database (str): name of OE database that 
        password (str): [description]
        user (str, optional): [description]. Defaults to 'root'.
        oe_url_path (str, optional): [description]. Defaults to 'http://localhost'.
    Raises:
        requests.exceptions.HTTPError: [description]
        RuntimeError: [description]
        RuntimeError: [description]
    """
    # Login info for OE database
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.13; rv:68.0) Gecko/20100101 Firefox/68.0'
    }
    login_data = {
        'db_name': database,
        'user': user,
        'password': password,
        'desired_action': 'login',
    }

    try:
        print('\nAccessing Open Enventory to generate structure images, please wait...!')

        # Get Login page and input username and password and login
        with requests.Session() as s:
            s.headers.update({'User-Agent': headers['User-Agent']})
            login_url = "{}/index.php".format(oe_url_path)

            # Login into the page
            # https://stackoverflow.com/a/21930636/6596203
            login = s.post(login_url, data=login_data)
            # soup = BeautifulSoup(login.content, "html.parser")
            # print(login.text)
            # print(soup.find_all('input'))

            # Throw error if cannot login
            if login.raise_for_status():
                print('Error: problems logging into OE!'.upper())
            else:
                # Call Open Enventory Batch processing to generate structure images
                data_to_be_fixed = {
                    'desired_action': 'fix_structures',
                    'save_settings': 'true',
                    'db_names[]': database,
                    'molecule': '1',
                    'emp_formula': '1',
                    'mw': '1',
                    'rdb': '1',
                    'smiles': '1',
                    'molfile_blob': '1',
                    'fingerprint': '1',
                }
                url2 = '{}/root_db_man.php'.format(oe_url_path)

                with Spinner("just waiting a bit.. \n"):
                    try:
                        second_request = s.post(url2, data=data_to_be_fixed)
                        # If the response was successful, no Exception will be raised
                        second_request.raise_for_status()
                    except requests.exceptions.HTTPError as http_err:
                        raise requests.exceptions.HTTPError(f'HTTP error occurred: {http_err}')  # Python 3.6
                    except Exception as err:
                        raise RuntimeError(f'Other error occurred: {err}')  # Python 3.6
                    else:
                        print('\nSuccess!')
                    time.sleep(1)

    except Exception as error:
        raise RuntimeError(error)


class Spinner:
    """To make spinning cursor while waiting:
    https://stackoverflow.com/a/58174909/6596203
    
    ### usage example:
    #
    # with Spinner("just waiting a bit.. "):
    #    do_something()
    #    time.sleep(3)

    """

    def __init__(self, message, delay=0.1):
        self.spinner = itertools.cycle(['-', '/', '|', '\\'])
        self.delay = delay
        self.busy = False
        self.spinner_visible = False
        sys.stdout.write(message)

    def write_next(self):
        with self._screen_lock:
            if not self.spinner_visible:
                sys.stdout.write(next(self.spinner))
                self.spinner_visible = True
                sys.stdout.flush()

    def remove_spinner(self, cleanup=False):
        with self._screen_lock:
            if self.spinner_visible:
                sys.stdout.write('\b')
                self.spinner_visible = False
                if cleanup:
                    sys.stdout.write(' ')       # overwrite spinner with blank
                    sys.stdout.write('\r')      # move to next line
                sys.stdout.flush()

    def spinner_task(self):
        while self.busy:
            self.write_next()
            time.sleep(self.delay)
            self.remove_spinner()

    def __enter__(self):
        if sys.stdout.isatty():
            self._screen_lock = threading.Lock()
            self.busy = True
            self.thread = threading.Thread(target=self.spinner_task)
            self.thread.start()

    def __exit__(self, exc_type, exc_val, exc_traceback):
        if sys.stdout.isatty():
            self.busy = False
            self.remove_spinner(cleanup=True)
        else:
            sys.stdout.write('\r')


def is_binary_string(bytes: bytes):
    '''To check if a file is a binary string
    Ref: https://stackoverflow.com/a/7392391/6596203
    Usage: 
        >>> is_binary_string(open('/usr/bin/python', 'rb').read(1024))
        True
        >>> is_binary_string(open('/usr/bin/dh_python3', 'rb').read(1024))
        False
    '''
    textchars = bytearray({7,8,9,10,12,13,27} | set(range(0x20, 0x100)) - {0x7f})
    return bool(bytes.translate(None, textchars))


def is_empty_mol_file(mol_file):
    '''To check if a downloaded mol file is empty.
    Some downloaded mol file are empty. Such file has content like this:
        '        
            SciTegic12151716442D

            0  0  0  0  0  0            999 V2000
            M  END
        '
    '''

    # '''Case 1: file contains just blank lines'''
    content = open(mol_file, 'r').read()
    # print(content)
    if re.search(r'^\s*$', content):
        return True

    '''Case 2: file contains blank mol file'''
    pattern = re.compile(r'^(?:\s*0)+\s*999 V2000$')
    with open(mol_file, 'r') as f:
        for line in f.readlines():
            result = pattern.search(line)
            if result:
                return True
        return False


if __name__ == '__main__':
    main()
