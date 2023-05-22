#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 10:49:17 2021

@author: earezza
"""

import os
import argparse
import pandas as pd
import numpy as np
import tqdm
import time
from openpyxl import load_workbook
from openpyxl.worksheet.worksheet import Worksheet

# DEFINE COMMANDLINE ARGUMENTS
describe_help = 'python organize_gene_results.py -f chromosome_1.xlsx -r result_filename.xlsx'
parser = argparse.ArgumentParser(description=describe_help)
parser.add_argument('-f', '--file', help='Full path to file for input', type=str)
parser.add_argument('-t', '--threshold', help='Threshold percentage for filtering', type=float, default=25)
parser.add_argument('-r', '--result', help='Full path to result file for output', type=str, default=os.getcwd() + '/organized_results.xlsx')
args = parser.parse_args()

# MAIN RUN
if __name__ == '__main__':
    t_start = time.time()
    # Show args for double-checking input
    print(args)
    
    # Read file
    print('Reading file...')
    df = pd.read_csv(args.file, na_values=['NO'])
    df.drop(index=df[df[df.columns[0]].isna()].index, inplace=True)
    #df.replace(to_replace=np.nan, value=0, inplace=True)
    
    num_cols = df.shape[1] - 1
    # Columns are either Gene name, PIPE soy-soy score
    # OR
    # Gene name, SPRINT, soy-soy score
    # OR
    # Gene name, PIPE soy-soy score, SPRINT soy-soy score
    # Columns also go in descending % (ranges 100% - 0%)
    
    print('Processing data...')
    # For file with ONLY PIPE or SPRINT columns
    if num_cols == 1:
        # Group genes based on desired criteria
        # Add gene for number of occurences in each criteria
        df_count = df.append(df[df[df.columns[1]] >= args.threshold])
        
        # Count number of occurrences
        counts = df_count.value_counts(subset=[df_count.columns[0]])
        counts = counts - 1
        gene_counts = counts.reset_index()
        
        # Passed threshold
        one = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 1][gene_counts.columns[0]])].sort_values(by=[df.columns[1]], ascending=False)
        # Didn't pass threshold
        zero = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 0][gene_counts.columns[0]])].sort_values(by=df.columns[1], ascending=False)

        groups = [zero, one]
        
    elif num_cols == 2:
        
        # Group genes based on desired criteria
        # Add gene for number of occurences in each criteria
        df_count = df.append(df[df[df.columns[1]] >= args.threshold])
        df_count = df_count.append(df[df[df.columns[1]] >= args.threshold])
        
        # Count number of occurrences
        counts = df_count.value_counts(subset=[df_count.columns[0]])
        counts = counts - 1
        gene_counts = counts.reset_index()
        
        # Passed both criteria
        two = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 2][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[2]], ascending=False)
        # PIPE soy-soy AND SPRINT soy-soy
        two = two[(two[two.columns[1]] >= args.threshold)
                        & (two[two.columns[2]] >= args.threshold)].sort_values(by=[df.columns[1], df.columns[2]], ascending=False)
        
        # Passed only one criteria
        one = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 1][gene_counts.columns[0]])].sort_values(by=[df.columns[1]], ascending=False)
        # Passed one (PIPE or SPRINT)
        one_first = one[(one[one.columns[1]] >= args.threshold) & (one[one.columns[2]] < args.threshold) ].sort_values(by=[df.columns[1], df.columns[2]], ascending=False)
        # Passed other (PIPE or SPRINT)
        one_second  = one[(one[one.columns[2]] >= args.threshold) & (one[one.columns[1]] < args.threshold) ].sort_values(by=[df.columns[2], df.columns[1]], ascending=False)
        # Combine
        one_all = pd.concat([one_first, one_second])
        
        # Didn't pass either criteria
        zero = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 0][gene_counts.columns[0]])].sort_values(by=df.columns[1], ascending=False)

        groups = [zero, one_all, two]

    else:
        print("Check number of columns in file")
        #exit()
    
    print('Creating organized results...')
    final = pd.DataFrame()
    for g in range(len(groups) - 1, -1, -1):
        final = final.append(pd.DataFrame(data={'Gene Name': ['PASSED %s CRITERIA (%s genes)'%(g, groups[g].shape[0])]}))
        final = final.append(groups[g])
        final = final.append(pd.DataFrame(data={'Gene Name': ['']}))
        
    
    # Output to file
    # Create file to append if not existing
    if not os.path.exists(args.result):
        writer = pd.ExcelWriter(args.result, engine='openpyxl')
        pd.DataFrame().to_excel(writer, sheet_name='Sheet1', index=False)
        writer.save()
        writer.close()
    
    print('Writing to file...')
    # Load file
    book = load_workbook(args.result)
    writer = pd.ExcelWriter(args.result, engine='openpyxl', mode='a')
    writer.book = book
    writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
    
    final.to_excel(writer, sheet_name='Sheet1', index=False)
    
    '''
    # Write grouping result
    for g in range(len(groups) - 1, -1, -1):
        groups[g].to_excel(writer, sheet_name='Sheet1', index=False)
        pd.DataFrame(columns=['Group %s'%g]).to_excel(writer, sheet_name='Sheet1', index=False)
    '''
    writer.save()
    writer.close()

    print('Done!')
    
    
    
