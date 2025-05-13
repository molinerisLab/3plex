#!/bin/env python
#import pandas as pd
#import numpy as np
#from scipy.sparse import dok_matrix
import argparse
from collections import defaultdict
import json
import msgpack
import sys
import os 
import numpy as np
import csv

def parse_csv(filename):
    quartiles = []; stability = []
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            quartiles.append(int(row[0]))
            stability.append(float(row[1]))
    return (quartiles, stability)

def get_p_value_for_tfo_regions(tpx_files):
    #1: Parse DBD regions => array con [Start, End]
    dbd_file = open(tpx_files[0],"r")
    dbds = []
    csvreader = csv.reader(dbd_file, delimiter='\t')
    for row in csvreader:
        dbds.append( (int(row[1]), int(row[2])) )
    dbd_file.close()
    
    #Read upper quartiles from stability file
    N, N_stab = parse_csv(tpx_files[1])
    #Read upper quartiles from all shuffled sequences
    all_dbds_formatted = []
    for random_file in tpx_files[2:]:
        all_dbds_formatted.append(parse_csv(random_file))

    
    #2: parsifichi i file con le randomizzazioni individualmente:
    #   =>Per ogni DBD, upper quartile di ogni randomizzazione
    #Contare in quale percentuale l'upper quartile è maggiore di N
    dbds_pvalue = []
    # [ ([a],[b])[ per ogni dbd]  per ogni random]
    for dbd_num, dbd in enumerate(dbds):
        dbd_formatted = [[list_[0][dbd_num],list_[1][dbd_num]] for list_ in all_dbds_formatted]
        #dbd_formatted = [get_upper_quartiles(random_tpx, [dbd]) for random_tpx in random_files]
        v = np.array([dbd_f[0] for dbd_f in dbd_formatted])
        count_N = sum(v>N[dbd_num]) / len(tpx_files[2:]) + (1/len(tpx_files[2:]))
        
        v = np.array([dbd_f[1] for dbd_f in dbd_formatted])
        count_N_stab = sum(v>N_stab[dbd_num]) / len(tpx_files[2:]) + (1/len(tpx_files[2:]))

        dbds_pvalue.append(f"{dbd[0]}\t{dbd[1]} \t{count_N}\t{count_N_stab}")
    return dbds_pvalue

def main():
    parser = argparse.ArgumentParser(description="read tpx stability file and compute tfo profiles")
    parser.add_argument('tpx_files', nargs='*', metavar='file', help='Input file(s)')
    args = parser.parse_args()

    result = get_p_value_for_tfo_regions(args.tpx_files)
    print("DBD_start\tDBD_end\ttts_count_pval\tstability_norm_pval")
    for r in result:
        print(r)
    
if __name__=="__main__":
    main()
    
