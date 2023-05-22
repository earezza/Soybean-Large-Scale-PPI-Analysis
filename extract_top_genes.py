#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 10:51:33 2021

Description:
    1. Read files in folder where each file contains one-to-all gene interactions and scores.
        i) Only keep longest-sequenced gene isoform if different isoforms exist in files (based on .fasta from --sequences input)
    2. For each file:
        i) Sort interactors in descending order by score (last column in each file)
        ii) Pull top X (20 for soy-scn, 40 for soy-soy) scoring interactors for given gene file (# top scorers is an option)
        iii) Calculate % of top scorers that include genes of interest
    3. Add top scorers to Excel under gene column
        i) Format and save as Excel, highlighting any genes of interest found in each column

Usage:
    Open terminal (command-line interface)
    Change to directory/folder where extract_top_genes.py is saved
    
    Run the following:
        python extract_top_genes.py -f <PATH_TO_FOLDER/> -r <path_to_result_filename> -t <number_of_top_scorers> -s <path_to_sequences_file>
    
        e.g. 
        python extract_top_genes.py -f Documents/SOY/soy-soy/ -r Documents/SOY/soy_top40 -t 40 -s Documents/SOY/Wm82.a2.v1.protein.aa
        
        Where soy-soy/ contains all gene .csv prediction files
        This will create an Excel with top 40 scorers saved as soy_top40.xlsx under Documents/SOY/
    
Requirements:
    This worked using the following (newer isoforms may also work):
    python 3.7.10
    pandas 0.24.2
    openpyxl 3.0.7
    tqdm 4.59.0

@author: Eric Arezza
Last updated: Nov. 3, 2021
"""

import os
import atexit
import argparse
import traceback
import pandas as pd
import numpy as np
import tqdm
import time
from openpyxl import load_workbook
from openpyxl.worksheet.worksheet import Worksheet

# DEFINE GENES OF INTEREST TO HIGHLIGHT RED IN EXCEL AND REPORT % FOUND IN TOP SCORERS
GENES_OF_INTEREST = [
    
    ]

MAX_SHEETS = 200
MAX_COLS = 800

# DEFINE COMMANDLINE ARGUMENTS
describe_help = 'python extract_top_genes.py -f PATH_TO_FOLDER/ -r path_to_result_filename.csv -t 40 -s sequences.fasta'
parser = argparse.ArgumentParser(description=describe_help)
parser.add_argument('-f', '--files', help='Full path to folder with .csv gene files', type=str)
parser.add_argument('-r', '--result', help='Full path to result file', type=str, default=os.getcwd() + '/top_genes.xlsx')
parser.add_argument('-t', '--top', help='Number of top interactors to include', type=int)
parser.add_argument('-s', '--sequences', help='Full path to fasta file containing all sequences', type=str)
parser.add_argument('-a', '--all', help='Flag to all gene isoforms', action='store_true')
parser.add_argument('-scn_only', '--scn_only', help='Flag to extract only scn scores from any soy gene files', action='store_true')
parser.add_argument('-soy_only', '--soy_only', help='Flag to extract only scn scores from any scn gene files', action='store_true')
args = parser.parse_args()

# DEFINE USEFUL FUNCTIONS
def get_isoforms(df_fasta):
    # Format df for easier use columns as ID, SEQUENCE
    df = df_fasta.copy()
    seq_id = df.iloc[::2, :].reset_index(drop=True)
    seq = df.iloc[1::2, :].reset_index(drop=True)
    seq_id.insert(1, 1, seq)
    seq_id[0] = seq_id[0].str.replace('>', '')
    seq_id[0] = seq_id[0].str.replace(r'.p', ' ', regex=True)
    seq_id[0] = [ i[0] for i in seq_id[0].str.split(' ') ]
    df = seq_id.copy()
    print('\t%s total genes in sequences file.'%df.shape[0])
    
    # Sort genes by having single or multi isoforms
    isoforms = np.array([ v[-1] for v in df[0].str.split('.') ])
    df.insert(2, 'Isoform', isoforms)
    df[0] = np.array([ '.'.join(i[:-1]) for i in df[0].str.split('.') ])
    seq_lengths = np.array([ len(i) for i in df[1] ])
    df = pd.DataFrame(data={0: df[0], 1: seq_lengths, 'Isoform': df['Isoform']})
    multi_isoforms = df[df.duplicated(subset=[0], keep=False)]
    multi_isoforms.reset_index(drop=True, inplace=True)
    single_isoforms = df[~df[0].isin(multi_isoforms[0].unique())]
    single_isoforms.reset_index(drop=True, inplace=True)
    
    # Sort multi isoforms by sequence length
    #multi_isoforms = multi_isoforms.sort_values(by=[0,1], ascending=False)
    multi_isoforms = multi_isoforms.sort_values(by=[0,'Isoform'], ascending=True)
    multi_isoforms.reset_index(inplace=True, drop=True)
    
    # Get multi isoforms with same length
    same_length = multi_isoforms[multi_isoforms.duplicated(subset=[0,1], keep=False)]
    same_length.reset_index(drop=True, inplace=True)
    
    # Get longest isoforms from genes with multiple isoforms
    longest = multi_isoforms.drop_duplicates(subset=[0], keep='first')
    
    # Get genes with multiple equally longest isoforms
    many_longest = same_length.merge(longest, on=[0, 1])
    many_longest.drop(columns='Isoform_y', inplace=True)
    many_longest.rename(columns={'Isoform_x': 'Isoform'}, inplace=True)
    many_longest.reset_index(drop=True, inplace=True)
    
    # Combine all longest isoforms and all genes with only one isoform, removing redundant
    isoforms = longest.append(many_longest, ignore_index=True)
    isoforms = isoforms.append(single_isoforms, ignore_index=True)
    isoforms = isoforms.drop_duplicates(subset=[0, 1])
    isoforms.reset_index(drop=True, inplace=True)
    isoforms.sort_values(by=[0, 'Isoform'], inplace=True)
    isoforms.reset_index(drop=True, inplace=True)
    isoforms[0] = isoforms[0] + '.' + isoforms['Isoform']
    
    many_longest[0] = many_longest[0] + '.' + many_longest['Isoform']

    return isoforms[0], many_longest[0]
    
def highlight_genes(gene):
    bg_color = 'pink' if pd.Series([gene]).isin(GENES_OF_INTEREST)[0] == True else 'none'
    return 'background-color: %s' % (bg_color)
def color_genes(gene):
    color = 'red' if pd.Series([gene]).isin(GENES_OF_INTEREST)[0] == True else 'black'
    return 'color: %s' % (color)

def get_saved_genes_sheet(filename, sheetname=None):
    # Start new Excel file if not found
    if sheetname == None:
        sheet = 'Sheet1'
    else:
        sheet = sheetname
    if not os.path.exists(filename):
        return np.array([], dtype=str), sheet
    
    # Get list of gene columns already recorded
    genes_recorded = np.array([], dtype=str)
    excel_file = pd.ExcelFile(filename)
    
    # Start new file if file sheets are maxed out
    if len(excel_file.sheet_names) > MAX_SHEETS:
        args.result = ''.join(args.result[:-1]) + '_new.' + args.result[-1]
    
    # Record saved gene columns
    for sheets in excel_file.sheet_names:
        genes_recorded = np.append(genes_recorded, excel_file.parse(sheets).columns)
    
    # Get current sheetname and number of columns used
    num_cols = excel_file.parse(sheet).shape[1]
    if num_cols < MAX_COLS:
        return genes_recorded, sheet
    else:
        return genes_recorded, "".join(filter(lambda x: not x.isdigit(), sheetname))+str(len(excel_file.sheet_names)+1)
    
def write_to_excel(file, df):
    saved_genes, sheet = get_saved_genes_sheet(file)
    # Create file to append if not existing
    if not os.path.exists(file):
        writer = pd.ExcelWriter(file, engine='openpyxl')
        pd.DataFrame().to_excel(writer, sheet_name=sheet, index=False)
        writer.save()
        writer.close()
        
    # Load file to append
    book = load_workbook(file)
    writer = pd.ExcelWriter(file, engine='openpyxl', mode='a')
    writer.book = book
    writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
    
    # Add SOY gene columns first
    #pbar = tqdm.tqdm(total=np.array(df[[ x for x in df.columns if 'Glyma' in x ]]).shape[0])
    pbar = tqdm.tqdm(total=df.columns.isin(saved_genes).shape[0] - df.columns.isin(saved_genes).sum())
    for i in range(0, 21):
        to_write, sheet = get_df_chromosome(df, num=i)
        if to_write.empty:
            continue
        excel = to_write.style.applymap(highlight_genes)
        excel = excel.applymap(color_genes)
        excel.to_excel(writer, sheet_name=sheet, index=False)
        saved_genes = np.append(saved_genes, to_write.columns)
        Worksheet(book, title=sheet)
        pbar.update(to_write.shape[1])
        df.drop(columns=to_write.columns, inplace=True)
        
    writer.save()
    writer.close()
    #pbar.close()
    
    # Load file to append
    book = load_workbook(file)
    writer = pd.ExcelWriter(file, engine='openpyxl', mode='a')
    writer.book = book
    writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
    
    # Add other gene columns until completed
    cols = 0
    sheet_num = 1
    sheet = 'Sheet1'
    #pbar = tqdm.tqdm(total=df.columns.isin(saved_genes).shape[0] - df.columns.isin(saved_genes).sum())
    while not all(df.columns.isin(saved_genes)):
        to_write = df[df.columns[cols:cols+MAX_COLS]]
        #to_write = get_df_chromosome(df, num=sheet_num)
        excel = to_write.style.applymap(highlight_genes)
        excel = excel.applymap(color_genes)
        excel.to_excel(writer, sheet_name=sheet, index=False)
        saved_genes = np.append(saved_genes, to_write.columns)
        sheet = "".join(filter(lambda x: not x.isdigit(), sheet))+str(sheet_num+1)
        Worksheet(book, title=sheet)
        pbar.update(to_write.shape[1])
        sheet_num += 1
        cols += MAX_COLS
        
    writer.save()
    writer.close()
    pbar.close()
    
def get_df_chromosome(df, num=0):
    if num == 0:
        chromosome_soy = df[[ c for c in df.columns if 'Glyma.U' in c ]]
        sheet = 'Soy Chromosome U'
    else:
        chromosome_soy = df[[ c for c in df.columns if ('Glyma.' +('%s'%num).zfill(2) + 'G') in c ]]
        sheet = 'Soy Chromosome %s'%num
    
    return chromosome_soy, sheet

def show_duration(t_start):
    # Display duration of run
    t_duration = time.time() - t_start
    day = t_duration // (24 * 3600)
    t_duration = t_duration % (24 * 3600)
    hour = t_duration // 3600
    t_duration %= 3600
    minutes = t_duration // 60
    t_duration %= 60
    seconds = t_duration
    print('Days', day)
    print('Hours', hour)
    print('Minutes', minutes)
    print('Seconds', seconds)

# MAIN RUN
if __name__ == '__main__':
    t_start = time.time()
    # Show args for double-checking input
    print(args)
    
    try:
    
        # If only want to get scn scores from files
        if args.scn_only:
            files = np.array([ f for f in os.listdir(path=args.files) if '.csv' in f and 'Glyma' in f ])
        # If only want to get soy scores from files
        elif args.soy_only:
            files = np.array([ f for f in os.listdir(path=args.files) if '.csv' in f and 'Heygly' in f ])
        else:
            # Get all gene files in folder
            files = np.array([ f for f in os.listdir(path=args.files) if '.csv' in f and ('Glyma' in f or 'Hetgly' in f) ])
        print('\nNumber of gene files:', len(files))
        
        # Skip any saved gene columns if exist already
        saved_genes, sheetname = get_saved_genes_sheet(args.result)
        print('%s genes already saved'%saved_genes.shape[0])
        
        if not args.all:
            print('Getting longest sequenced isoforms...')
            seq = pd.read_csv(args.sequences, sep='\n', header=None)
            isoforms, many_longest = get_isoforms(seq)
            print('\t%s relevant gene isoforms\n\t%s have multiple equally long sequences...'%(isoforms.shape[0], isoforms[isoforms.isin(many_longest)].shape[0]))
            # Only process gene files for relevant isoforms
            #files = np.array([ i for i in files if i.replace('.csv', '' ) in isoforms.values or i.replace('.csv', '' ) in many_longest.values ])
            files = np.array([ i for i in files if i.replace('.csv', '' ) in isoforms.values and i.replace('.csv', '' ) not in saved_genes ])
            print('\n%s relevant files to process...\n'%len(files))
        
        # For result file
        final = pd.DataFrame()
        final_isoforms = pd.DataFrame()
        
        # Iterate through each file
        print('\nIterating through files...')
        for f in tqdm.tqdm(files):
    
            # Get gene name from filename
            gene = f.replace('.csv', '')
            
            # Read file
            df = pd.read_csv(args.files + f)
            
            # Remove any genes that are not SOY or SCN
            df = df[(df[df.columns[1]].str.contains('Glyma')) | (df[df.columns[1]].str.contains('Hetgly'))]
            df.reset_index(drop=True, inplace=True)
            
            # If only want to consider SCN for interactor scores
            if args.scn_only:
                df = df[df[df.columns[1]].str.contains('Hetgly')]
                df.reset_index(drop=True, inplace=True)
            if args.soy_only:
                df = df[df[df.columns[1]].str.contains('Glyma')]
                df.reset_index(drop=True, inplace=True)
            
            # Sort by scores in descending order
            df.sort_values(by=df.columns[-1], ascending=False, inplace=True)
            df.reset_index(drop=True, inplace=True)
            
            # Keep copy for recording interactor isoform number
            df_isoforms = df.copy()
            # Remove gene isoform number (after the last '.')
            df[df.columns[1]] = [ '.'.join(g[:-1]) for g in df[df.columns[1]].str.split('.').values ]
            
            # Drop any duplicated interactor genes from list (removes isoform variants, keeping the top scored one)
            df = df.drop_duplicates(subset=[df.columns[1]])
            
            # Re-sort non-duplicated genes
            df.sort_values(by=df.columns[-1], ascending=False, inplace=True)
            
            # Get top interactors
            df = df.iloc[:args.top]
            df_isoforms = df_isoforms.iloc[df.index]
            
            df.reset_index(drop=True, inplace=True)
            df_isoforms.reset_index(drop=True, inplace=True)

            # Count % of top interactors found in genes of interest
            percent_interested = ( df[df.columns[1]].isin(GENES_OF_INTEREST).sum() / args.top )*100
            
            # Create column for gene and include % of top in genes of interest
            gene_info = pd.DataFrame(data=df[df.columns[1]].values, columns=[gene])
            gene_info_isoforms = pd.DataFrame(data=df_isoforms[df_isoforms.columns[1]].values, columns=[gene])
            
            gene_info = gene_info[gene].append(pd.Series(percent_interested), ignore_index=True)
            gene_info_isoforms = gene_info_isoforms[gene].append(pd.Series(percent_interested), ignore_index=True)
            
            gene_info = pd.DataFrame(gene_info, columns=[gene])
            gene_info_isoforms = pd.DataFrame(gene_info_isoforms, columns=[gene])
    
            #Add gene info to final
            final.insert(len(final.columns), gene, gene_info[gene])
            final_isoforms.insert(len(final_isoforms.columns), gene, gene_info_isoforms[gene])
        
        final = final[sorted(final.columns)]
        final_isoforms = final_isoforms[sorted(final_isoforms.columns)]
        show_duration(t_start)
        
        if final.empty:
            print('No columns to add...\nDone!\n')
            show_duration(t_start)
            exit()
        
        # Write to excel
        print('Creating coloured Excel, %s gene columns...'%final.shape[1])
        write_to_excel(args.result, final)
        
        show_duration(t_start)
        
        print('Creating Excel with isoform numbers, %s gene columns...'%final_isoforms.shape[1])
        iso_filename = args.result.split('.')
        iso_filename = ''.join(iso_filename[:-1]) + '_isoform_numbers.' + iso_filename[-1]
        write_to_excel(iso_filename, final_isoforms)
        
        print('Done!\n')
        show_duration(t_start)
    except (KeyboardInterrupt, Exception):
        atexit.register(show_duration, t_start)
        traceback.print_exc()
    
