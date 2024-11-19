#!/bin/env python
import numpy as np
from collections import defaultdict
import argparse
import msgpack
import sys
import csv

def get_upper_quartiles(tpx_file, dbds):
    profile = defaultdict(lambda: 0)
    profile_stab = defaultdict(lambda: [])
    FILE_LEN = sum(1 for line in tpx_file)
    tpx_file.seek(0)
    for index, line in enumerate(tpx_file):
        b, e, stability = line.rstrip().split("\t")
        b=int(b)
        e=int(e)
        stability=float(stability)
        for i in range(b,e):
            profile[i] = profile[i]+1
            profile_stab[i] = profile_stab[i] + [stability]
        if index%500==0:
            print(f"P1: {index}/{FILE_LEN}", file=sys.stderr)
    quartiles = []
    quartiles_stability = []
    for index, dbd in enumerate(dbds):
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
        if index%100==0:
            print(f"DBD: {index}/{len(dbds)}", file=sys.stderr)
    return quartiles, quartiles_stability

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="read tpx stability file and compute tfo profiles")
    parser.add_argument('tpx_files', nargs='*', metavar='file', help='Input file(s)')
    args = parser.parse_args()

    dbd_file = open(args.tpx_files[0],"r")
    dbds = []
    csvreader = csv.reader(dbd_file, delimiter='\t')
    for row in csvreader:
        dbds.append( (int(row[1]), int(row[2])) )
    dbd_file.close()
    with open(args.tpx_files[1], "r") as tpx_file:
        result = get_upper_quartiles(tpx_file, dbds)
    sys.stdout.buffer.write(msgpack.packb(result))