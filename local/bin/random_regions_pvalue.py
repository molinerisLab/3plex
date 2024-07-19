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

def get_p_value_for_tfo_regions(tpx_files):
    def get_upper_quartiles(tpx_file, dbds):
        profile = defaultdict(lambda: 0)
        profile_stab = defaultdict(lambda: [])
        for line in tpx_file:
            b, e, stability = line.rstrip().split("\t")
            b=int(b)
            e=int(e)
            stability=float(stability)
            for i in range(b,e):
                profile[i] = profile[i]+1
                profile_stab[i] = profile_stab[i] + [stability]
        quartiles = []
        quartiles_stability = []
        for dbd in dbds:
            values = []
            values_stability = []
            for i in range(dbd[0], dbd[1]):
                values.append(profile[i])
                values_stability += profile_stab[i]
            
            upper_quartile = np.percentile(values, 0.75)
            quartiles.append(upper_quartile)
            
            if len(values_stability) > 0:
                upper_quartile_stability = np.percentile(values_stability, 0.75)
            else:
                upper_quartile_stability = 0
            quartiles_stability.append(upper_quartile_stability)
        tpx_file.seek(0)
        return quartiles, quartiles_stability
        
    dbd_file = open(tpx_files[0],"r")
    stability_file = open(tpx_files[1],"r")
    random_files = [open(f,"r") for f in tpx_files[2:]]

    #1: Parse DBD regions => array con [Start, End]
    dbds = []
    csvreader = csv.reader(dbd_file, delimiter='\t')
    for row in csvreader:
        dbds.append( (int(row[1]), int(row[2])) )
    #2: For each DBD, upper quartile in stability_file
    N, N_stab = get_upper_quartiles(stability_file, dbds)
    #2: parsifichi i file con le randomizzazioni individualmente:
    #   =>Per ogni DBD, upper quartile di ogni randomizzazione
    #Contare in quale percentuale l'upper quartile è maggiore di N
    dbds_pvalue = []
    all_dbds_formatted = [get_upper_quartiles(random_tpx, dbds) for random_tpx in random_files]
    # [ ([a],[b])[ per ogni dbd]  per ogni random]
    for dbd_num, dbd in enumerate(dbds):
        dbd_formatted = [[list_[0][dbd_num],list_[1][dbd_num]] for list_ in all_dbds_formatted]
        #dbd_formatted = [get_upper_quartiles(random_tpx, [dbd]) for random_tpx in random_files]
        v = [dbd_f[0] for dbd_f in dbd_formatted]
        count_N = sum(v>N[dbd_num]) / len(random_files) + (1/len(random_files))
        
        v = [dbd_f[1] for dbd_f in dbd_formatted]
        count_N_stab = sum(v>N_stab[dbd_num]) / len(random_files) + (1/len(random_files))

        dbds_pvalue.append(f"{dbd[0]}\t{dbd[1]} \t{count_N}\t{count_N_stab}")
    dbd_file.close()
    stability_file.close()
    for r in random_files:
        r.close()
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
    
