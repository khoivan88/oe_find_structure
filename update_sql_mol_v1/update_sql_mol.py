#!/usr/bin/python
import mysql.connector as mariadb

def main():

    mariadb_connection, cursor = connect('root', 'Buchwa1dr0ck5', 'test1')
    search_empty_smile(cursor, 'limit 10')
    mariadb_connection.close()

def connect(username, pw, database):
    # Create a cursor:
    # https://mariadb.com/resources/blog/how-connect-python-programs-mariadb
    mariadb_connection = mariadb.connect(user=username, password=pw, database=database)
    cursor = mariadb_connection.cursor()
    return mariadb_connection, cursor

def search_empty_smile(cursor, limit=''):
    cursor.execute("SELECT cas_nr, molecule_id, smiles, molfile_blob FROM molecule where smiles='' {}".format(limit))

    for row in cursor:
        print('CAS# {}: {}'.format(row[0], row))

if __name__ == '__main__':
    main()
