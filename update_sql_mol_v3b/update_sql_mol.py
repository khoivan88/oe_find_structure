'''
This version 3b: different than version 3. This version download all of the CAS
files in the molecule table. The reason is because OE when updating structures
after the mol files are uploaded onto mysql get rid of all of the existing structures.

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
'''


#!/usr/bin/python
import mysql.connector as mariadb

import requests
import wget
from pathlib import Path
from multiprocessing import Pool
# from functools import partial
import getpass

missing_mol_file = set()
download_path = '/var/lib/mysql/missing_mol_files/'

def main():
    global download_path

    # Get user input for root password and the database needs to be updated
    # to hide password input: https://stackoverflow.com/questions/9202224/getting-command-line-password-input-in-python
    # password = input('Please type in the password for "root" user: ')
    password = getpass.getpass('Please type in the password for "root" user: ')
    database = input('Please type in the name of the database needs updating: ')
    # mariadb_connection, cursor = connect('root', 'Buchwa1dr0ck5', 'test1')

    # Info for mysql connection and query can be found here:
    # https://mariadb.com/resources/blog/how-connect-python-programs-mariadb
    # https://dev.mysql.com/doc/connector-python/en/connector-python-tutorial-cursorbuffered.html
    mariadb_connection = mariadb.connect(user='root', password=password, database=database)

    cursor_select = mariadb_connection.cursor(buffered=True)

    print('Getting molecule with missing structures. Please wait!')
    query = ("SELECT distinct cas_nr FROM molecule WHERE cas_nr!=''")

    try:
        cursor_select.execute(query)
    except mariadb.Error as error:
        print('Error: {}'.format(error))


    to_be_downloaded = []
    for (cas_nr, ) in cursor_select:
        to_be_downloaded.append(cas_nr)
    # print('There are {} molecules missing structure'.format(len(to_be_downloaded)))

    print('Downloading missing mol files. Please wait!')
    try:
        with Pool(100) as p:
            p.map(extracting_mol, to_be_downloaded)
    except Exception as error:
        print(error)

    finally:
        print('Updating SQL table!')
        for cas in to_be_downloaded:
            try:
                update_sql(mariadb_connection, cas)
            except mariadb.Error as error:
                print('Error: {}'.format(error))

    mariadb_connection.close()
    print('\n{} mol files are missing: '.format(len(missing_mol_file)))
    print(missing_mol_file)

def update_sql(mariadb_connection, cas_nr):
    global download_path
    cursor_update = mariadb_connection.cursor(buffered=True)
    mol_file = Path(download_path + '{}.mol'.format(cas_nr))
    # print(mol_file)

    # result = extracting_mol(cas_nr)
    # if molfile exists or downloaded (extracting_mol return -1 or 0)
    if mol_file.exists():
        print('CAS# {}'.format(cas_nr))
        # cursor_update.execute(insert_mol_file, (mol_file, cas_nr))
        cursor_update.execute("UPDATE molecule SET molfile_blob=LOAD_FILE('{}') WHERE cas_nr='{}'".format(mol_file, cas_nr))
        mariadb_connection.commit()
        print('\tmol file uploaded successfully!')
    # extracting_mol return the cas# of those that it could not find mol file
    else:
        missing_mol_file.add(cas_nr)


def extracting_mol(cas_nr):
    global download_path
    '''This function is used to extract a single mol file
    See here for more info: http://stackabuse.com/download-files-with-python/'''


    # print('Beginning download file with wget')

    headers = {
        'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/53.0.2785.143 Safari/537.36'}

    #get url from ChemicalBook to download mol file
    url = 'http://www.chemicalbook.com/CAS/mol/'
    cas = cas_nr
    file_name = cas + '.mol'
    full_url = url + file_name

    # download_path = '/Users/khoivan/Downloads/mol_files/'
    download_file = Path(download_path + file_name)

    # Check if the file not exists and download
    #check file exists: https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists
    if download_file.exists():
        # print('{} already downloaded'.format(file_name))
        return -1
    else:
        # check if the mol find exist. requests.get().history would show code
        # 302 (redirect) if the mol file does not exist.
        # ChemicalBook automatically redirect if could not find mol file
        if (len(requests.get(full_url, headers=headers, timeout=5).history) > 0):
            #return the cas # of the mol whose mol file cannot be found
            return cas
        else:
            print('Downloading {} ...'.format(file_name))
            wget.download(full_url, download_path[:-1])
            return 0


def connect(username, pw, database):
    # Create a cursor:
    # https://mariadb.com/resources/blog/how-connect-python-programs-mariadb
    return mariadb_connection, cursor

if __name__ == '__main__':
    main()
