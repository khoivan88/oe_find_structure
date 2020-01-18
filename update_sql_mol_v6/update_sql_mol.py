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
"""

# from functools import partial
import getpass
import itertools
import os
import sys
import threading
import time
from multiprocessing import Pool
from pathlib import Path

import mysql.connector as mariadb
import requests
import wget

# from bs4 import BeautifulSoup


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

        # Step1: run SELECT query to find CAS#
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

        # Step 2: downloading mol file
        # Check if download path with the missing_mol_file directory exists. If not, create it
        # https://stackoverflow.com/questions/12517451/automatically-creating-directories-with-file-output
        # https://docs.python.org/3/library/os.html#os.makedirs
        os.makedirs(download_path, exist_ok=True)

        print('Downloading missing mol files. Please wait!')
        try:
            with Pool(25) as p:
                p.map(extracting_mol, to_be_downloaded)
        except Exception as error:
            print(error)

        # Step 3: run UPDATE query to upload
        finally:
            # Upload mol files into mySQL
            print('Updating SQL table!')
            count_file_updated = 0
            # for (cas, ) in cursor_select:
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

            # Login into OE and use Batch Processing to generate structure images
            oe_url_path = 'http://localhost' if oe_url_path_input == '' else oe_url_path_input
            oe_batch_process(
                database=database,
                user='root',
                password=password,
                oe_url_path=oe_url_path)

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

    # result = extracting_mol(cas_nr)
    # if molfile exists or downloaded (extracting_mol return -1 or 0)
    if mol_file.exists():
        print('CAS# {:<15}: '.format(cas_nr), end='')
        # cursor_update.execute(insert_mol_file, (file_path, cas_nr))
        cursor_update.execute("UPDATE molecule SET molfile_blob=LOAD_FILE('{}') WHERE cas_nr='{}'".format(mol_file, cas_nr))
        mariadb_connection.commit()
        print('mol file uploaded successfully!')
        return 1
    # extracting_mol return the cas# of those that it could not find mol file
    else:
        missing_mol_file.add(cas_nr)
        return 0


def extracting_mol(cas_nr):
    global download_path
    '''
    This function is used to extract a single mol file
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
    download_file = Path(download_path + '/' + file_name)

    # Check if the file not exists and download
    # check file exists: https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists
    if download_file.exists():
        # print('{} already downloaded'.format(file_name))
        return -1
    else:
        # check if the mol find exist. requests.get().history would show code
        # 302 (redirect) if the mol file does not exist.
        # ChemicalBook automatically redirect if could not find mol file
        # redirected = False;
        try:
            hist_len = len(requests.get(full_url, headers=headers, timeout=10).history)
            if hist_len > 0:  # requests.get.history will a list of more than 0 element if redirected
                # redirected = True
                # return the cas # of the mol whose mol file cannot be found
                return cas
            else:
                print('Downloading {} ...'.format(file_name))
                wget.download(full_url, download_path, bar=wget.bar_thermometer)
                return 0

        except Exception as error:
            print(error)
            # return 1


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


if __name__ == '__main__':
    main()
