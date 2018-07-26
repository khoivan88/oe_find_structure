#!/usr/bin/python
import mysql.connector as mariadb

def main():
    mariadb_connection = mariadb.connect(user='root', password='Buchwa1dr0ck5', database='test1')
    cursor = mariadb_connection.cursor()
    

if __name__ == '__main__':
    main()
