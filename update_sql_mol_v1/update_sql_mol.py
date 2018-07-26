#!/usr/bin/python
import mysql.connector as mariadb

def main():
    mariadb_connection = mariadb.connect(user='root', password='Buchwa1dr0ck5', database='test1')
    cursor = mariadb_connection.cursor()
    cursor.execute("SELECT cas_nr, molecule_id, smiles, molfile_blob FROM molecule limit 50")

    for row in cursor:
        print('CAS# {}: {}'.format(row[0], row))

if __name__ == '__main__':
    main()
