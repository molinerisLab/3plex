#!/bin/env python
import argparse
from collections import defaultdict
import json
import msgpack
import sys
import os 
import numpy as np


# stdin has 3 columns:
#
#  1       TFO_start
#  2       TFO_end
#  3       Stability
#
# without header
#
# Output is a dict with 2 list
# respectively, "all" contains statistics computed on all tpx (puts together tpx from all random runs)
# while "best" takes, for each position, the best stability of each random run and computes statistics on it
# Both outputs are list, for each position of ssRNA they contains {avg, std}

def get_stability_profile(tpx_file):
    # assume data sorted on desending Stability

    #profiles=dok_matrix((stability_max, ssRNA_len),dtype=np.uintc)
    profile = defaultdict(lambda: [])
    
    for line in tpx_file:
        b, e, stability = line.rstrip().split("\t")
        b=int(b)
        e=int(e)
        stability=float(stability)
        for i in range(b,e):
            profile[i] = profile[i] + [stability]
    return profile


def profiles_multiple_stats(profiles_multiple):
    # Create an empty dictionary to store the cumulative sum of profile heights
    profile_matrix_all = defaultdict(lambda: [])

    # Iterate over each profile in profiles_multiple
    for profile in profiles_multiple:
        for position in profile.keys():
            profile_matrix_all[position] = profile_matrix_all[position] + profile[position]
    #profile_matrix_all contains the list of stability of all tpxs on all randomizations
    #Now compute avg, stdev
    statistics_all = {}
    for position in profile_matrix_all.keys():
        statistics_all[position] = {
            "avg": np.mean(profile_matrix_all[position]),
            "std": np.std(profile_matrix_all[position])
        }
    #Repeat but with best stability for each random run
    profile_matrix_best = defaultdict(lambda: [])
    # Iterate over each profile in profiles_multiple
    for profile in profiles_multiple:
        for position in profile.keys():
            profile_matrix_best[position] = profile_matrix_all[position] + [max(profile[position])]
    #profile_matrix_all contains the list of best stability for each random run
    #Now compute avg, stdev
    statistics_best = {}
    for position in profile_matrix_best.keys():
        statistics_best[position] = {
            "avg": np.mean(profile_matrix_best[position]),
            "std": np.std(profile_matrix_best[position])
        }
    return statistics_all, statistics_best

def dict_to_list_of_positions(profile):
    position_max = max(map(lambda x: int(x), list(profile.keys())))
    positions = [{"avg":0, "std":0} for _ in range(position_max+1)]
    for x in profile.keys():
        positions[int(x)] = profile[x]
    return positions

def main():
    parser = argparse.ArgumentParser(description="read tpx stability file and compute tfo profiles")
    parser.add_argument('tpx_files', nargs='*', metavar='file', help='Input file(s)')
    args = parser.parse_args()
    profiles_multiple=[]
    for file_name in args.tpx_files:
        with open(file_name) as my_file:
            data = get_stability_profile(my_file)
        profiles_multiple.append(data)
    stats_all, stats_best = profiles_multiple_stats(profiles_multiple)
    to_export = {
        "all": dict_to_list_of_positions(stats_all),
        "best": dict_to_list_of_positions(stats_best)
    }
    sys.stdout.write(json.dumps(to_export))    

if __name__=="__main__":
    main()