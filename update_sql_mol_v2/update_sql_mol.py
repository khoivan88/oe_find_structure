#!/usr/bin/python
import mysql.connector as mariadb
from pathlib import Path

import requests
import wget
from pathlib import Path
from multiprocessing import Pool
from functools import partial

missing_mol_file = set()
download_path = '/var/lib/mysql/missing_mol_files/'

def main():
    global download_path
    # mariadb_connection, cursor = connect('root', 'Buchwa1dr0ck5', 'test1')
    # search_empty_smile(cursor, 'limit 1')
    # search_empty_smile(cursor)
    mariadb_connection = mariadb.connect(user='root', password='Buchwa1dr0ck5', database='test1')

    cursor_select = mariadb_connection.cursor(buffered=True)

    # query = ("SELECT cas_nr FROM molecule WHERE smiles='' limit 1")
    query = ("SELECT distinct cas_nr FROM molecule WHERE smiles=''")
    cursor_select.execute(query)

    to_be_downloaded = []
    for (cas_nr, ) in cursor_select:
        to_be_downloaded.append(cas_nr)
    # print('There are {} molecules missing structure'.format(len(to_be_downloaded)))

    try:
        with Pool(10) as p:
            p.map(extracting_mol, to_be_downloaded)
            # missing_mol_file.append([x for x in record if x is not None])

    finally:
    # func = partial(update_sql, mariadb_connection)
        for cas in to_be_downloaded:
            update_sql(mariadb_connection, cas)

    mariadb_connection.close()
    print('\n{} mol files are missing: '.format(len(missing_mol_file)))
    print(missing_mol_file)

def update_sql(mariadb_connection, cas_nr):
    global download_path
    cursor_update = mariadb_connection.cursor(buffered=True)
    print('CAS# {}'.format(cas_nr))
    mol_file = Path(download_path + '{}.mol'.format(cas_nr))
    # print(mol_file)

    # result = extracting_mol(cas_nr)
    # if molfile exists or downloaded (extracting_mol return -1 or 0)
    if mol_file.exists():
        # cursor_update.execute(insert_mol_file, (mol_file, cas_nr))
        cursor_update.execute("UPDATE molecule SET molfile_blob=LOAD_FILE('{}') WHERE cas_nr='{}'".format(mol_file, cas_nr))
        mariadb_connection.commit()
        print('\n\tmol file uploaded successfully!')
    # extracting_mol return the cas# of those that it could not find mol file
    else:
        missing_mol_file.add(cas_nr)


def some_func():
    try:
        #Reading csv file: see more https://automatetheboringstuff.com/chapter14/
        exampleFile = open(file)
        exampleReader = csv.reader(exampleFile)
        exampleData = list(exampleReader)

        # from: https://medium.com/python-pandemonium/how-to-speed-up-your-python-web-scraper-by-using-multiprocessing-f2f4ef838686
        # Pool(10) means 10 entries will be processed at a single time
        with Pool(10) as p:
            record = p.map(extracting_mol, exampleData)
            missing_mol_file.append([x for x in record if x is not None])

        # if record is not None (None means it found the url) then add to missing mol_file list:

        print('\n\nCould NOT find mol file for: {}'.format(missing_mol_file))

    except Exception as ex:
        print(ex)

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

# def search_empty_smile(cursor, limit=''):
#     # cursor.execute("SELECT cas_nr, molecule_id, smiles, molfile_blob FROM molecule where smiles='' {}".format(limit))
#     cursor.execute("SELECT cas_nr FROM molecule where smiles='' {}".format(limit))
#     # missing_structure_cas_list = cursor[:]
#
#     for row in cursor:
#         print('CAS# {}'.format(row))
#
#     print('cursor type is: {}'.format(type(cursor)))
#     # print(len(missing_structure_cas_list))

def insert_mol(cursor, limit=''):
    try:
        cursor.execute("SELECT cas_nr FROM molecule where smiles='' {}".format(limit))

        cas_nr = cursor[0][0]

        print (cas_nr)
        # print('cursor type is: {}'.format(type(cursor)))
        # print(len(missing_structure_cas_list))
    except mariadb.Error as error:
        print('Error: {}'.format(error))

if __name__ == '__main__':
    main()
