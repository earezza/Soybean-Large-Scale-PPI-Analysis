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
parser.add_argument('-r', '--result', help='Full path to result file for output', type=str, default=os.getcwd() + '/organized_results.xlsx')
args = parser.parse_args()

'''
# When PIPE SOY-SCN is not considered/passed and < 25, sort by SPRINT SOY-SCN instead
def resort(df):
    above = df[df[df.columns[4]] >= 25]
    below = df[df[df.columns[4]] < 25 ]
    below = below.sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False)
    df_resorted = above.append(below)
    return df_resorted
'''

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
    # Genes with Yes's for SNP column go first, then LOF column, then PIPE soy-SCN, then annotation, then PIPE soy-soy
    # Columns also go in descending % for each part for the soy-SCN column as well
    # range for the soy-SCN column is anything from 100-25% is considered a yes. obviously the ones closest to 100% should be at the top as seen in the screenshot
    # range for the soy-soy column is anything from 100-30% is considered a yes
    
    print('Processing data...')
    # For file with ONLY PIPE or SPRINT columns
    if num_cols == 5:
        # Group genes based on desired criteria
        # Gene Name, SNPs, LOFs, PIPE Soy-Soy, PIPE Soy-SCN, Annotations
        # Add gene for occurence in each criteria
        df_count = df.append(df[~df[df.columns[1]].isna()])
        df_count = df_count.append(df[~df[df.columns[2]].isna()])
        df_count = df_count.append(df[df[df.columns[4]] >= 25])
        df_count = df_count.append(df[~df[df.columns[5]].isna()])
        df_count = df_count.append(df[df[df.columns[3]] >= 30])
        
        # Count number of occurrences
        counts = df_count.value_counts(subset=[df_count.columns[0]])
        counts = counts - 1
        gene_counts = counts.reset_index()
        
        # Group by number of criteria passed and sorted descending based on priority of columns (as discussed)
        five = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 5][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[4], df.columns[3], df.columns[2], df.columns[5]], ascending=False).replace(to_replace=np.nan, value='NO')
       
        four = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 4][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[4], df.columns[3], df.columns[2], df.columns[5]], ascending=False)
        # Passed all except LOF
        four_first = four[four[four.columns[2]].isna()].sort_values(by=[df.columns[4], df.columns[3], df.columns[1], df.columns[5]], ascending=False).sort_values(by=[df.columns[4], df.columns[3], df.columns[1] , df.columns[5]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed all except SNPs
        four_second = four[four[four.columns[1]].isna()].sort_values(by=[df.columns[4], df.columns[3], df.columns[1], df.columns[5]], ascending=False).sort_values(by=[df.columns[4], df.columns[3], df.columns[2], df.columns[5]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed all except SOY-SOY
        four_third = four[four[four.columns[3]] < 30].sort_values(by=[df.columns[4], df.columns[3], df.columns[1], df.columns[5]], ascending=False).sort_values(by=[df.columns[4], df.columns[1], df.columns[5], df.columns[2]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed all except SOY-SCN
        four_fourth = four[four[four.columns[4]] < 25].sort_values(by=[df.columns[4], df.columns[3], df.columns[1], df.columns[5]], ascending=False).sort_values(by=[df.columns[3], df.columns[1], df.columns[2], df.columns[5]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed all except annotation
        four_fifth = four[four[four.columns[5]].isna()].sort_values(by=[df.columns[4], df.columns[3], df.columns[1], df.columns[5]], ascending=False).sort_values(by=[df.columns[4], df.columns[3], df.columns[1], df.columns[2]], ascending=False).replace(to_replace=np.nan, value='NO')
        four_all = pd.concat([four_first, four_second, four_third, four_fourth, four_fifth])

        three = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 3][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[4], df.columns[3], df.columns[2], df.columns[5]], ascending=False)
        # Passed SNPs, SOY-SCN, Annotation
        three_first = three[(~three[three.columns[1]].isna()) 
                         & (~three[three.columns[5]].isna()) 
                         & (three[three.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[3], df.columns[1]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed LOF, SOY-SCN, Annotation
        three_second = three[(~three[three.columns[2]].isna()) 
                         & (~three[three.columns[5]].isna()) 
                         & (three[three.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[3], df.columns[2]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SNPs, LOF, Annotation
        three_third = three[(~three[three.columns[1]].isna()) 
                         & (~three[three.columns[2]].isna()) 
                         & (~three[three.columns[5]].isna())].sort_values(by=[df.columns[4], df.columns[3], df.columns[1], df.columns[2]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed LOF, SOY-SOY, Annotation
        three_fourth = three[(~three[three.columns[2]].isna()) 
                         & (~three[three.columns[5]].isna()) 
                         & (three[three.columns[3]] >= 30)].sort_values(by=[df.columns[4], df.columns[3], df.columns[1], df.columns[2]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SOY-SCN, SOY-SOY, Annotation
        three_fifth = three[(~three[three.columns[5]].isna()) 
                         & (three[three.columns[4]] >= 25)
                         & (three[three.columns[3]] >= 30)].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SNPs and SOY-SCN and LOFs
        three_sixth = three[(~three[three.columns[1]].isna()) 
                        & (~three[three.columns[2]].isna()) 
                        & (three[three.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[3], df.columns[1]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed LOF and SOY-SOY and SOY-SCN
        three_seventh = three[(~three[three.columns[2]].isna()) 
                              & (three[three.columns[4]] >= 25)
                              & (three[three.columns[3]] >= 30)].sort_values(by=[df.columns[4], df.columns[3], df.columns[1]], ascending=False).replace(to_replace=np.nan, value='NO')     
        # Passed SNPs and SOY-SOY and annotations
        three_eighth = three[(~three[three.columns[1]].isna()) 
                             & (~three[three.columns[5]].isna()) 
                             & (three[three.columns[3]] >= 30)].sort_values(by=[df.columns[4], df.columns[3], df.columns[1]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SNPs, LOFs, SOY-SOY
        three_ninth = three[(~three[three.columns[1]].isna()) 
                            & (~three[three.columns[2]].isna()) 
                            & (three[three.columns[3]] >= 30)].sort_values(by=[df.columns[4], df.columns[3], df.columns[1]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SNPs, SOY-SOY, SOY-SCN
        three_tenth = three[(~three[three.columns[1]].isna()) 
                            & (three[three.columns[3]] >= 30) 
                            & (three[three.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[3], df.columns[1]], ascending=False).replace(to_replace=np.nan, value='NO')
        three_all = pd.concat([three_first, three_second, three_third, three_fourth, three_fifth, three_sixth, three_seventh, three_eighth, three_ninth, three_tenth])
        #leftover_three = three[~three[three.columns[0]].isin(three_all[three_all.columns[0]])].sort_values(by=[df.columns[4], df.columns[3], df.columns[1], df.columns[2], df.columns[5]], ascending=False).replace(to_replace=np.nan, value='NO')
        #three_all = three_all.append(leftover_three)
        #three_all.drop_duplicates(subset=[three_all.columns[0]], inplace=True)
        
        two = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 2][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[4], df.columns[3], df.columns[2], df.columns[5]], ascending=False)
        # Passed SNPs, Annotation
        two_first = two[ (~two[two.columns[1]].isna())
                         & (~two[two.columns[5]].isna())].sort_values(by=[df.columns[4], df.columns[3], df.columns[1]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed LOF, Annotation
        two_second = two[(~two[two.columns[2]].isna()) 
                         & (~two[two.columns[5]].isna())].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SNPs, SOY-SCN
        two_third = two[(~two[two.columns[1]].isna()) 
                         & (two[two.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed LOF, SOY-SCN
        two_fourth = two[(~two[two.columns[2]].isna())
                         & (two[two.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SOY-SCN, Annotation
        two_fifth = two[(two[two.columns[4]] >= 25) 
                         & (~two[two.columns[5]].isna())].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SOY-SOY, Annotation
        two_sixth = two[(two[two.columns[3]] >= 30) 
                         & (~two[two.columns[5]].isna())].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SNPs and SOY-SOY
        two_seventh = two[(~two[two.columns[1]].isna())
                         & (two[two.columns[3]] >= 30)].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed LOF and SOY-SOY
        two_eighth = two[(~two[two.columns[2]].isna())
                         & (two[two.columns[3]] >= 30)].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SNPs and LOFs
        two_ninth = two[(~two[two.columns[1]].isna())
                         & (~two[two.columns[2]].isna())].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SOY-SOY and SOY-SCN
        two_tenth = two[(two[two.columns[3]] >= 30)
                         & (two[two.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')

        two_all = pd.concat([two_first, two_second, two_third, two_fourth, two_fifth, two_sixth, two_seventh, two_eighth, two_ninth, two_tenth])
        #leftover_two = two[~two[two.columns[0]].isin(two_all[two_all.columns[0]])].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        #two_all = two_all.append(leftover_two)
        #two_all.drop_duplicates(subset=[two_all.columns[0]], inplace=True)
        
        one = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 1][gene_counts.columns[0]])].sort_values(by=[df.columns[3], df.columns[4], df.columns[5], df.columns[1], df.columns[2]], ascending=False)
        # Passed SNPs
        one_first = one[(~one[one.columns[1]].isna())].sort_values(by=[df.columns[1], df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed LOF
        one_second = one[(~one[one.columns[2]].isna())].sort_values(by=[df.columns[4], df.columns[3], df.columns[2]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed Annotation
        one_third = one[(~one[one.columns[5]].isna())].sort_values(by=[df.columns[4], df.columns[3], df.columns[5]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SOY-SCN
        one_fourth = one[(one[one.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Passed SOY-SOY
        one_fifth = one[(one[one.columns[3]] >= 30)].sort_values(by=[df.columns[3], df.columns[4]], ascending=False).replace(to_replace=np.nan, value='NO')
        one_all = pd.concat([one_first, one_second, one_third, one_fourth, one_fifth])

        zero = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 0][gene_counts.columns[0]])].sort_values(by=[df.columns[4], df.columns[3]], ascending=False).replace(to_replace=np.nan, value='NO')

        groups = [zero, one_all, two_all, three_all, four_all, five]
    
    # For file with PIPE AND SPRINT columns
    elif num_cols == 7:
        # Gene Name, SNPs, LOFs, PIPE Soy-Soy, PIPE Soy-SCN, SPRINT Soy-Soy, SPRINT Soy-SCN, Annotations
        # Add gene for occurence in each criteria
        df_count = df.append(df[~df[df.columns[1]].isna()])
        df_count = df_count.append(df[~df[df.columns[2]].isna()])
        df_count = df_count.append(df[df[df.columns[4]] >= 25])
        df_count = df_count.append(df[~df[df.columns[7]].isna()])
        df_count = df_count.append(df[df[df.columns[3]] >= 30])
        df_count = df_count.append(df[df[df.columns[5]] >= 30])
        df_count = df_count.append(df[df[df.columns[6]] >= 25])
        
        # Count number of occurrences
        counts = df_count.value_counts(subset=[df_count.columns[0]])
        counts = counts - 1
        gene_counts = counts.reset_index()
        
        
        # Group by number of criteria passed and sorted descending based on priority of columns (as discussed)
        # Sort by PIPE soy-scn
        seven = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 7][gene_counts.columns[0]])].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        
        six = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 6][gene_counts.columns[0]])].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False)
        # All except LOF
        six_first = six[(six[six.columns[2]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SNP
        six_second = six[(six[six.columns[1]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SPRINT soy-soy
        six_third = six[six[six.columns[5]] < 30].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SPRINT soy-scn
        six_fourth = six[six[six.columns[6]] < 25].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except PIPE soy-soy
        six_fifth = six[six[six.columns[3]] < 30].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except PIPE soy-scn
        six_sixth = six[six[six.columns[4]] < 25].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except annotations
        six_seventh = six[six[six.columns[7]].isna()].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Combine all six
        six_all = pd.concat([six_first, six_second, six_third, six_fourth, six_fifth, six_sixth, six_seventh])
        
        
        five = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 5][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[4], df.columns[6], df.columns[3], df.columns[5], df.columns[2], df.columns[5]], ascending=False)
        # All except LOF and SPRINT soy-soy
        five_first = five[(five[five.columns[2]].isna())
                        & (five[five.columns[5]] < 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except LOF and PIPE soy-soy
        five_second = five[(five[five.columns[2]].isna())
                        & (five[five.columns[3]] < 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SNP and SPRINT soy-soy
        five_third = five[(five[five.columns[1]].isna())
                        & (five[five.columns[5]] < 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SNP and PIPE soy-soy
        five_fourth = five[(five[five.columns[1]].isna())
                        & (five[five.columns[3]] < 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SPRINT soy-scn and SPRINT soy-soy
        five_fifth = five[(five[five.columns[5]] < 30)
                        & (five[five.columns[6]] < 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except PIPE soy-scn and PIPE soy-soy
        five_sixth = five[(five[five.columns[3]] < 30)
                        & (five[five.columns[4]] < 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SNP and SPRINT soy-scn
        five_seventh = five[(five[five.columns[1]].isna())
                        & (five[five.columns[6]] < 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SNP and PIPE soy-scn
        five_eighth = five[(five[five.columns[1]].isna())
                        & (five[five.columns[4]] < 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SNP and LOF
        five_ninth = five[(five[five.columns[1]].isna())
                        & (five[five.columns[2]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SPRINT soy-soy and annotations
        five_tenth = five[(five[five.columns[7]].isna())
                        & (five[five.columns[5]] < 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except PIPE soy-soy and annotations
        five_eleventh = five[(five[five.columns[7]].isna())
                        & (five[five.columns[3]] < 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except LOF annotations
        five_twelfth = five[(five[five.columns[2]].isna())
                        & (five[five.columns[7]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SNP and annotations
        five_thirteenth = five[(five[five.columns[1]].isna())
                        & (five[five.columns[7]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except LOF and SPRINT soy-scn
        five_fourteenth = five[(five[five.columns[2]].isna())
                        & (five[five.columns[6]] < 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except LOF and PIPE soy-scn
        five_fifteenth = five[(five[five.columns[2]].isna())
                        & (five[five.columns[4]] < 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except SPRINT soy-scn and annotations
        five_sixteenth = five[(five[five.columns[7]].isna())
                        & (five[five.columns[6]] < 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except PIPE soy-scn and annotations
        five_seventeenth = five[(five[five.columns[7]].isna())
                        & (five[five.columns[4]] < 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except PIPE soy-soy and SPRINT soy-soy
        five_eighteenth = five[(five[five.columns[3]] < 30)
                        & (five[five.columns[5]] < 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except PIPE soy-soy and SPRINT soy-scn
        five_nineteenth = five[(five[five.columns[3]] < 30)
                        & (five[five.columns[6]] < 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except PIPE soy-scn and SPRINT soy-soy
        five_twenty = five[(five[five.columns[4]] < 25)
                        & (five[five.columns[5]] < 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # All except PIPE soy-scn and SPRINT soy-scn
        five_twentyone = five[(five[five.columns[4]] < 25)
                        & (five[five.columns[6]] < 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Combine all five
        five_all = pd.concat([five_first, five_second, five_third, five_fourth, five_fifth, five_sixth, five_seventh, five_eighth,
                              five_ninth, five_tenth, five_eleventh, five_twelfth, five_thirteenth, five_fourteenth,
                              five_fifteenth, five_sixteenth, five_seventeenth, five_eighteenth, five_nineteenth,
                              five_twenty, five_twentyone])
        
        
        four = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 4][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[4], df.columns[6], df.columns[3], df.columns[5], df.columns[2], df.columns[5]], ascending=False)
        # SNPs, PIPE soy-soy, PIPE soy-scn, annotations
        four_first = four[(~four[four.columns[1]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, SPRINT soy-soy, SPRINT soy-scn, annotations
        four_second = four[(~four[four.columns[1]].isna())
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, PIPE soy-scn, annotations
        four_third = four[(~four[four.columns[2]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, SPRINT soy-soy, SPRINT soy-scn, annotations
        four_fourth = four[(~four[four.columns[2]].isna())
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, PIPE soy-scn, annotations
        four_fifth = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, SPRINT soy-scn, annotations
        four_sixth = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[6]] >= 25)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, PIPE soy-soy, PIPE soy-scn
        four_seventh = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, SPRINT soy-soy, SPRINT soy-scn
        four_eighth = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-scn, SPRINT soy-scn, annotations
        four_ninth = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[7]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy, SPRINT soy-scn, annotations
        four_tenth = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[7]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy, SPRINT soy-soy, annotations
        four_eleventh = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[7]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-scn, SPRINT soy-soy, annotations
        four_twelfth = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[7]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-scn, SPRINT soy-scn, annotations
        four_thirteenth = four[(~four[four.columns[2]].isna())
                        & (~four[four.columns[7]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, SPRINT soy-scn, annotations
        four_fourteenth = four[(~four[four.columns[2]].isna())
                        & (~four[four.columns[7]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, SPRINT soy-soy, annotations
        four_fifteenth = four[(~four[four.columns[2]].isna())
                        & (~four[four.columns[7]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-scn, SPRINT soy-soy, annotations
        four_sixteenth = four[(~four[four.columns[2]].isna())
                        & (~four[four.columns[7]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, PIPE soy-scn, SPRINT soy-scn, annotations
        four_seventeenth = four[(~four[four.columns[7]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-scn, SPRINT soy-soy, SPRINT soy-scn, annotations
        four_eighteenth = four[(~four[four.columns[7]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-scn, SPRINT soy-soy, SPRINT soy-scn
        four_nineteenth = four[(~four[four.columns[1]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy, SPRINT soy-soy, SPRINT soy-scn
        four_twenty = four[(~four[four.columns[1]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy, PIPE soy-scn, SPRINT soy-scn
        four_twentyone = four[(~four[four.columns[1]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy, PIPE soy-scn, SPRINT soy-soy
        four_twentytwo = four[(~four[four.columns[1]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-scn, SPRINT soy-soy, SPRINT soy-scn
        four_twentythree = four[(~four[four.columns[2]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, SPRINT soy-soy, SPRINT soy-scn
        four_twentyfour = four[(~four[four.columns[2]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, PIPE soy-scn, SPRINT soy-scn
        four_twentyfive = four[(~four[four.columns[2]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, PIPE soy-scn, SPRINT soy-soy
        four_twentysix = four[(~four[four.columns[2]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, PIPE soy-scn, SPRINT soy-soy, SPRINT soy-scn
        four_twentyseven = four[(four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, PIPE soy-scn, SPRINT soy-scn
        four_twentyeight = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, PIPE soy-soy, SPRINT soy-scn
        four_twentynine = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, PIPE soy-scn, SPRINT soy-soy
        four_thirty = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, PIPE soy-soy, annotations
        four_thirtyone = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, PIPE soy-soy, SPRINT soy-soy
        four_thirtytwo = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[3]] >= 30)
                        & (four[four.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, SPRINT soy-soy, annotations
        four_thirtythree = four[(~four[four.columns[1]].isna())
                        & (~four[four.columns[2]].isna())
                        & (four[four.columns[5]] >= 30)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, SPRINT soy-soy, SPRINT soy-scn, annotations
        four_thirtyfour = four[(four[four.columns[3]] >= 30)
                        & (four[four.columns[5]] >= 30)
                        & (four[four.columns[6]] >= 25)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, PIPE soy-scn, SPRINT soy-soy, annotations
        four_thirtyfive = four[(four[four.columns[3]] >= 30)
                        & (four[four.columns[4]] >= 25)
                        & (four[four.columns[5]] >= 30)
                        & (~four[four.columns[7]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Combine all four
        four_all = pd.concat([four_first, four_second, four_third, four_fourth, four_fifth, four_sixth, four_seventh, four_eighth,
                              four_ninth, four_tenth, four_eleventh, four_twelfth, four_thirteenth, four_fourteenth,
                              four_fifteenth, four_sixteenth, four_seventeenth, four_eighteenth, four_nineteenth,
                              four_twenty, four_twentyone, four_twentytwo, four_twentythree, four_twentyfour,
                              four_twentyfive, four_twentysix, four_twentyseven, four_twentyeight, four_twentynine,
                              four_thirty, four_thirtyone, four_thirtytwo, four_thirtythree, four_thirtyfour,
                              four_thirtyfive])
        
        
        three = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 3][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[4], df.columns[6], df.columns[3], df.columns[5], df.columns[2], df.columns[5]], ascending=False)
        # SNPs, PIPE soy-scn, annotations
        three_first = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[4]] >= 25)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, SPRINT soy-scn, annotations
        three_second = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[6]] >= 25)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-scn, annotations
        three_third = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[4]] >= 25)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, SPRINT soy-scn, annotations
        three_fourth = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[6]] >= 25)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, annotations
        three_fifth = three[(~three[three.columns[1]].isna())
                        & (~three[three.columns[2]].isna())
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, annotations
        three_sixth = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, SPRINT soy-soy, annotations
        three_seventh = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[5]] >= 30)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, PIPE soy-scn, annotations
        three_eighth = three[(three[three.columns[3]] >= 30)
                        & (three[three.columns[4]] >= 25)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SPRINT soy-soy, SPRINT soy-scn, annotations
        three_ninth = three[(three[three.columns[5]] >= 30)
                        & (three[three.columns[6]] >= 25)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, PIPE soy-scn
        three_tenth = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[4]] >= 25)
                        & (~three[three.columns[2]].isna())].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, SPRINT soy-scn
        three_eleventh = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[6]] >= 25)
                        & (~three[three.columns[2]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy, annotations
        three_twelfth = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, SPRINT soy-soy, annotations
        three_thirteenth = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[5]] >= 30)
                        & (~three[three.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy, PIPE soy-scn
        three_fourteenth = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (three[three.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, SPRINT soy-soy, SPRINT soy-scn
        three_fifteenth = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[5]] >= 30)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, PIPE soy-scn
        three_sixteenth = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (three[three.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, SPRINT soy-soy, SPRINT soy-scn
        three_seventeenth = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[5]] >= 30)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, PIPE soy-soy
        three_eighteenth = three[(~three[three.columns[1]].isna())
                        & (~three[three.columns[2]].isna())
                        & (three[three.columns[3]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs, SPRINT soy-soy
        three_nineteenth = three[(~three[three.columns[1]].isna())
                        & (~three[three.columns[2]].isna())
                        & (three[three.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-scn, SPRINT soy-scn, annotations
        three_twenty = three[(~three[three.columns[7]].isna())
                        & (three[three.columns[4]] >= 25)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, SPRINT soy-soy, annotations (this was 8th # PIPE soy-soy, PIPE soy-scn, annotations)
        three_twentyone = three[(~three[three.columns[7]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (three[three.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-scn, SPRINT soy-soy, annotations
        three_twentytwo = three[(~three[three.columns[7]].isna())
                        & (three[three.columns[4]] >= 25)
                        & (three[three.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, SPRINT soy-scn, annotations (this was 9th # SPRINT soy-soy, SPRINT soy-scn, annotations)
        three_twentythree = three[(~three[three.columns[7]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-scn, SPRINT soy-scn
        three_twentyfour = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[4]] >= 25)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-scn, SPRINT soy-soy
        three_twentyfive = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[4]] >= 25)
                        & (three[three.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy, SPRINT soy-scn
        three_twentysix = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy, SPRINT soy-soy
        three_twentyseven = three[(~three[three.columns[1]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (three[three.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-scn, SPRINT soy-scn
        three_twentyeight = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[4]] >= 25)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-scn, SPRINT soy-soy
        three_twentynine = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[4]] >= 25)
                        & (three[three.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, SPRINT soy-scn
        three_thirty = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy, SPRINT soy-soy
        three_thirtyone = three[(~three[three.columns[2]].isna())
                        & (three[three.columns[3]] >= 30)
                        & (three[three.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, PIPE soy-scn, SPRINT soy-scn
        three_thirtytwo = three[(three[three.columns[3]] >= 30)
                        & (three[three.columns[4]] >= 25)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, PIPE soy-scn, SPRINT soy-soy
        three_thirtythree = three[(three[three.columns[3]] >= 30)
                        & (three[three.columns[4]] >= 25)
                        & (three[three.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, SPRINT soy-soy, SPRINT soy-scn
        three_thirtyfour = three[(three[three.columns[3]] >= 30)
                        & (three[three.columns[5]] >= 30)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-scn, SPRINT soy-soy, SPRINT soy-scn
        three_thirtyfive = three[(three[three.columns[4]] >= 25)
                        & (three[three.columns[5]] >= 30)
                        & (three[three.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Combine all three
        three_all = pd.concat([three_first, three_second, three_third, three_fourth, three_fifth, three_sixth, three_seventh,
                             three_eighth, three_ninth, three_tenth, three_eleventh, three_twelfth, three_thirteenth,
                             three_fourteenth, three_fifteenth, three_sixteenth, three_seventeenth,
                             three_eighteenth, three_nineteenth, three_twenty, three_twentyone, three_twentytwo,
                             three_twentythree, three_twentyfour, three_twentyfive, three_twentysix, three_twentyseven,
                             three_twentyeight, three_twentynine, three_thirty, three_thirtyone, three_thirtytwo,
                             three_thirtythree, three_thirtyfour, three_thirtyfive])
        
        two = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 2][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[4], df.columns[6], df.columns[3], df.columns[5], df.columns[2], df.columns[5]], ascending=False)
        # SNPs, and annotations
        two_first = two[(~two[two.columns[1]].isna())
                        & (~two[two.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, annotations
        two_second = two[(~two[two.columns[2]].isna())
                        & (~two[two.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-scn
        two_third = two[(~two[two.columns[1]].isna())
                        & (two[two.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, SPRINT soy-scn
        two_fourth = two[(~two[two.columns[1]].isna())
                        & (two[two.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-scn
        two_fifth = two[(~two[two.columns[2]].isna())
                        & (two[two.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, SPRINT soy-scn
        two_sixth = two[(~two[two.columns[2]].isna())
                        & (two[two.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-scn, annotations
        two_seventh = two[(~two[two.columns[7]].isna())
                        & (two[two.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SPRINT soy-scn, annotations
        two_eighth = two[(~two[two.columns[7]].isna())
                        & (two[two.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, annotations
        two_ninth = two[(~two[two.columns[7]].isna())
                        & (two[two.columns[3]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SPRINT soy-soy, annotations
        two_tenth = two[(~two[two.columns[7]].isna())
                        & (two[two.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, PIPE soy-soy
        two_eleventh = two[(~two[two.columns[1]].isna())
                        & (two[two.columns[3]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, SPRINT soy-soy
        two_twelfth = two[(~two[two.columns[1]].isna())
                        & (two[two.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, PIPE soy-soy
        two_thirteenth = two[(~two[two.columns[2]].isna())
                        & (two[two.columns[3]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs, SPRINT soy-soy
        two_fourteenth = two[(~two[two.columns[2]].isna())
                        & (two[two.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-scn, PIPE soy-soy
        two_fifteenth = two[(two[two.columns[3]] >= 30)
                        & (two[two.columns[4]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SPRINT soy-scn, SPRINT soy-soy
        two_sixteenth = two[(two[two.columns[5]] >= 30)
                        & (two[two.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SNPs, LOFs
        two_seventeenth = two[(~two[two.columns[1]].isna())
                        & (~two[two.columns[2]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-scn, SPRINT soy-scn
        two_eighteenth = two[(two[two.columns[4]] >= 25)
                        & (two[two.columns[6]] >= 25)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-scn, SPRINT soy-soy
        two_nineteenth = two[(two[two.columns[4]] >= 25)
                        & (two[two.columns[5]] >= 30)].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, SPRINT soy-scn
        two_twenty = two[(two[two.columns[3]] >= 30)
                        & (two[two.columns[6]] >= 25)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy, SPRINT soy-soy
        two_twentyone = two[(two[two.columns[3]] >= 30)
                        & (two[two.columns[5]] >= 30)].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Combine all two
        two_all = pd.concat([two_first, two_second, two_third, two_fourth, two_fifth, two_sixth, two_seventh,
                             two_eighth, two_ninth, two_tenth, two_eleventh, two_twelfth, two_thirteenth,
                             two_fourteenth, two_fifteenth, two_sixteenth, two_seventeenth, two_eighteenth,
                             two_nineteenth, two_twenty, two_twentyone])
        
        one = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 1][gene_counts.columns[0]])].sort_values(by=[df.columns[1], df.columns[4], df.columns[6], df.columns[3], df.columns[5], df.columns[2], df.columns[5]], ascending=False)
        # SNPs
        one_first = one[(~one[one.columns[1]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # LOFs
        one_second =  one[(~one[one.columns[2]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Annotations
        one_third = one[(~one[one.columns[7]].isna())].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-scn
        one_fourth = one[one[one.columns[4]] >= 25].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SPRINT soy-scn
        one_fifth = one[one[one.columns[6]] >= 25].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # PIPE soy-soy
        one_sixth = one[one[one.columns[3]] >= 30].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # SPRINT soy-soy
        one_seventh = one[one[one.columns[5]] >= 30].sort_values(by=[df.columns[6], df.columns[4], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        # Combine all one
        one_all = pd.concat([one_first, one_second, one_third, one_fourth, one_fifth, one_sixth, one_seventh])
        
        zero = df[df[df.columns[0]].isin(gene_counts[gene_counts[gene_counts.columns[-1]] == 0][gene_counts.columns[0]])].sort_values(by=[df.columns[4], df.columns[6], df.columns[1], df.columns[3], df.columns[5], df.columns[2], df.columns[7]], ascending=False).replace(to_replace=np.nan, value='NO')
        
        groups = [zero, one_all, two_all, three_all, four_all, five_all, six_all, seven]
        
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
    
    
    