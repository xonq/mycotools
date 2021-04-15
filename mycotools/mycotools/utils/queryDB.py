#! /usr/bin/env python3

import pandas as pd
import argparse, sys
from mycotools.lib.dbtools import query, db_edit, db2df, df2db, masterDB, formatPath



def main():

    parser = argparse.ArgumentParser(description='Query/edit columns.')
    parser.add_argument('-c', '--column', default = None, help='Only edit this column')
    parser.add_argument('-o','--ome',default='internal_ome', help='Index column. DEFAULT: internal_ome')
    parser.add_argument('-d','--database', default = masterDB(), help='Kontools db. DEFAULT: master')

    args = parser.parse_args()

    exit_code = 0
    print( '\n' + args.database , flush = True)
    if formatPath(args.database) == formatPath(masterDB()):
        print( 'ATTENTION: Querying master database' , flush = True)
    db_df = db2df(args.database)

    while exit_code == 0:
        info = input('\n\t[S]ave, [q]uery, [e]dit: ')
        if info in { 'q', 'Q' }:
            while exit_code == 0:
                query_list = []
                query = None
                while query != '':
                    query = input('Input ome or blank to query/exit: ')
                    query_list.append(query)
                print('\n', flush = True)
                if args.column in db_df.columns.values:
                    for query in query_list:
                        if query != '':
                            exit_code = query(db_df, 'ome', query, query2 = args.column, ome_index = args.ome)
                else:
                    for query in query_list:
                        if query != '':
                            exit_code = query(db_df, 'ome', query, ome_index = args.ome)
                if query_list == ['']:
                    break 
        elif info in { 'e', 'E' }:
            while exit_code == 0:
                query = input('Input ome or blank to edit/exit: ')
                if query == '':
                    break
                if args.column not in db_df.columns.values:
                    query2 = input('\nInput column to edit: ')
                    query(db_df, 'ome', query, query2 = query2, ome_index = args.ome)
                    edit = input('\nInput edit - blank to exit, space for blank edit: ')
                    if edit == '':
                        break
                    elif edit == ' ':
                        edit = ''
                    db_df, exit_code = db_edit(query, query2, edit, db_df, ome_index = args.ome)
                else:
                    query(db_df, 'ome', query, args.column, ome_index = args.ome)
                    edit = input('\nInput edit - blank to exit, space for blank edit: ')  
                    if edit == '':
                        break
                    elif edit == ' ':
                        edit = ''
                    db_df, exit_code = db_edit(query, args.column, edit, db_df, ome_index = args.ome)
        elif info in {'S', 's', ''}:
            exit_code = 1
 
    df2db( db_df, args.database )

if __name__ == '__main__':
    main()
    sys.exit(0)
