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

"""


#!/usr/bin/python

import mysql.connector as mariadb

import requests
import os
import wget
from pathlib import Path
from multiprocessing import Pool
# from functools import partial
import getpass

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
        query = ("SELECT distinct cas_nr FROM molecule WHERE cas_nr!='' and molfile_blob like '%open enventory%'")
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
            with Pool(100) as p:
                p.map(extracting_mol, to_be_downloaded)
        except Exception as error:
            print(error)

        # Step 3: run UPDATE query to upload
        finally:
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
            print(missing_mol_file)
            print('\nSummary: ')
            print('\t{} mol files are missing: '.format(len(missing_mol_file)))
            # path, dirs, files = next(os.walk(download_path))
            # file_count = len(files)
            print('\t{} mol files updated! '.format(count_file_updated))

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
    mol_file = Path(download_path + '/{}.mol'.format(cas_nr))
    # print(mol_file)

    # result = extracting_mol(cas_nr)
    # if molfile exists or downloaded (extracting_mol return -1 or 0)
    if mol_file.exists():
        print('CAS# {}'.format(cas_nr), end='')
        # cursor_update.execute(insert_mol_file, (mol_file, cas_nr))
        cursor_update.execute("UPDATE molecule SET molfile_blob=LOAD_FILE('{}') WHERE cas_nr='{}'".format(mol_file, cas_nr))
        mariadb_connection.commit()
        print('\tmol file uploaded successfully!')
        return 1
    # extracting_mol return the cas# of those that it could not find mol file
    else:
        missing_mol_file.add(cas_nr)
        return 0


def extracting_mol(cas_nr):
    global download_path
    '''This function is used to extract a single mol file
    See here for more info: http://stackabuse.com/download-files-with-python/'''

    headers = {
        'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/53.0.2785.143 Safari/537.36'}

    #get url from ChemicalBook to download mol file
    url = 'http://www.chemicalbook.com/CAS/mol/'
    cas = cas_nr
    file_name = cas + '.mol'
    full_url = url + file_name


    # download_path = '/Users/khoivan/Downloads/mol_files/'
    download_file = Path(download_path + '/' + file_name)

    # Check if the file not exists and download
    #check file exists: https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists
    if download_file.exists():
        # print('{} already downloaded'.format(file_name))
        return -1
    else:
        # check if the mol find exist. requests.get().history would show code
        # 302 (redirect) if the mol file does not exist.
        # ChemicalBook automatically redirect if could not find mol file
        redirected = False;
        try:
            hist_len = len(requests.get(full_url, headers=headers, timeout=20).history)
            if hist_len > 0:
                redirected = True
        except Exception as error:
            print(error)
            # return 1
        finally:
            if not redirected:
                print('Downloading {} ...'.format(file_name))
                wget.download(full_url, download_path, bar=wget.bar_thermometer)
                return 0
            else:
                #return the cas # of the mol whose mol file cannot be found
                return cas


# def connect(username, pw, database):
#     # Create a cursor:
#     # https://mariadb.com/resources/blog/how-connect-python-programs-mariadb
#     return mariadb_connection, cursor

if __name__ == '__main__':
    main()
