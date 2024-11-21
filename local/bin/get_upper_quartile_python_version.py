#!/bin/env python
import numpy as np
from collections import defaultdict
import argparse
import msgpack
import sys
import csv
import time


"""
    This is NOT the get_upper_quartile script used in the Snakefile rules.
    The right script is the one written in C. This one I leave as a "backup" but it's 4x slower
"""


def find_dbds(b, e, DBDS):
    to_return = []
    L, H = 0, len(DBDS)
    while(L <= H):
        M = (L+H)//2
        dbd = DBDS[M]
        if ((dbd[2] <= b and dbd[3] > b) or (dbd[2] < e and dbd[2] >= b)):
            to_return.append(dbd)
            #Check if other dbds just before or after ours are included
            i = M-1
            while i >= 0 and ((DBDS[i][2] <= b and DBDS[i][3] > b) or (DBDS[i][2] < e and DBDS[i][2] >= b)):
                to_return.append(DBDS[i])
                i -= 1 
            i = M+1
            while i < len(DBDS) and ((DBDS[i][2] <= b and DBDS[i][3] > b) or (DBDS[i][2] < e and DBDS[i][2] >= b)):
                to_return.append(DBDS[i])
                i += 1
            return to_return
        elif dbd[2] >= e: 
            H = M-1
        else:
            L = M+1
    return to_return

def get_upper_quartiles(tpx_file, dbds):
    #upper quartile computed on n. tpx for each position of the dbd
    # Only way is keeping count for each position of dbd.
    #upper quartile stability computed on all stability values of all tpx of the dbd
    # No need to track every position of the dbd, stability values can be kept in a single array,
    # but, they must be repeated for the number of positions, or the information must be tracked in some way
    #Es: (Value, Count). Sort by value, iterate once to get Sum(Count),iterate second time to get quartile :)
    
    #Benchmarks
    #TIME = time.time()
    TIME_FINDING = 0
    TIME_ITERATING = 0
    TIME_ITERATING_2 = 0

    #Each dbd as: [ [pos1, pos2, pos3, ...], [(stability1, count_1),  ...], begin, end ]
    DBDS = [[ [0 for _ in range(dbd[0], dbd[1])], [], dbd[0], dbd[1] ] for dbd in dbds]
    for line in tpx_file:
        b, e, stability = line.rstrip().split("\t")
        b=int(b)
        e=int(e)
        stability=float(stability)
        #Find DBDS where the interval is included
        TIME = time.time()
        current_dbds = find_dbds(b, e, DBDS)
        TIME_FINDING += time.time() - TIME
        TIME = time.time()
        for dbd in current_dbds:
            b_dbd = max(b - dbd[2],0) 
            e_dbd = min(max(e - dbd[2],0), dbd[3]-dbd[2])
            for i in range(b_dbd, e_dbd):
                dbd[0][i] += 1
            dbd[1].append((stability, e_dbd-b_dbd))
        TIME_ITERATING += time.time() - TIME
    
    #Compute quartiles
    quartiles = []
    quartiles_stability = []
    TIME = time.time()
    for dbd in DBDS:
        values = dbd[0]
        quartiles.append(np.percentile(values, 75))

        values_stability = sorted(dbd[1], key=lambda x: x[0])
        if len(values_stability) == 0:
            quartiles_stability.append(0)
            continue
        
        count = sum([v[1] for v in values_stability])
        threshold = count * 0.25
        for s in values_stability:
            count -= s[1]
            if count <= threshold:
                quartiles_stability.append(s[0])
                break
    TIME_ITERATING_2 = time.time() - TIME
    #sys.stderr.write(f"Finding: {TIME_FINDING}\nIterating: {TIME_ITERATING} - It2 {TIME_ITERATING_2}\n")
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
    for i in range(len(result[0])):
        print(f"{result[0][i]};{result[1][i]}")
    #sys.stdout.buffer.write(msgpack.packb(result))