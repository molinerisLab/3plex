#!/bin/env python
#import pandas as pd
#import numpy as np
#from scipy.sparse import dok_matrix
import argparse
from collections import defaultdict
import json
import sys
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
# assume data sorted on desending Stability
#
# output is a dictionary, keys are stability levels, values are the profiles of tfo counts


def get_TFO_profile_allSparse(tpx_file):
    # assume data sorted on desending Stability

    #profiles=dok_matrix((stability_max, ssRNA_len),dtype=np.uintc)
    profiles = {}
    profile_current = defaultdict(lambda: 0)
    stability_pre=None
    best_stability = {}
    max_len = 0
    
    for line in tpx_file:
        b, e, stability = line.rstrip().split("\t")
        b=int(b)
        e=int(e)
        if stability != stability_pre:
            if stability_pre is not None:
                profiles[stability_pre]=dict(profile_current)
            stability_pre=stability

        for i in range(b,e):
            profile_current[i]+=1
            if (i not in best_stability):
                best_stability[i] = stability
            if (e > max_len):
                max_len = e
                
    return {"profiles": profiles, "best_stability": best_stability}, max_len

class Ranges:
    container=[]
    def reset_container(self):
        self.container = []

    def add_range(self,b,e,count, i=0):
        if e<b:
            raise ValueError("invalid range %d-%d, count: %d, i: %d" % (b,e,count,i))
        #print("range %d-%d, count: %d, i: %d" % (b,e,count,i))
        if (b!=e):
            l = [count, (b + ((e-b)/2)), e-b+1]
        else:
            l = [count, b]
        self.container.append(l)

def profile2ranges(profiles):    
    profiles_range={}
    
    for stability in profiles.keys():
        profile = profiles[stability]
        profile_range = Ranges()
        profile_range.reset_container()
        
        range_b=None
        range_count=None
        i_pre = None
        keys = list(profile.keys())
        keys.sort(key=int)
        for i in keys:
            count = profile[i]
            i = int(i);
            count = int(count)
            if range_b is None:
                range_b = i
                i_pre = i
                range_count=count

            if count != range_count or i > i_pre + 1:
                profile_range.add_range(range_b, i_pre, range_count)
                range_b = i
                range_count = count
            i_pre = i
        
        profile_range.add_range(range_b,i_pre,range_count)
        #print("%d\t%d\t%d" % (range_b,i_pre,count))
        profiles_range[stability]= profile_range.container

    return(profiles_range)            

def best_stability_to_range(best_stability):
    profile_range = Ranges()
    profile_range.reset_container()
        
    range_b=None
    range_value=None
    i_pre = None
    keys = list(best_stability.keys())
    keys.sort(key=int)
    for i in keys:
        value = best_stability[i]
        i = int(i);
        if range_b is None:
            range_b = i
            i_pre = i
            range_value = value

        if value != range_value or i > i_pre + 1:
            profile_range.add_range(range_b, i_pre, range_value)
            range_b = i
            range_value = value
        i_pre = i
        
        profile_range.add_range(range_b,i_pre,range_value)
    return profile_range.container

def best_stability_to_array(best_stability, length):
    array = np.zeros(length)
    for key in best_stability.keys():
        #key to integer
        index = int(key)
        array[index] = best_stability[key]
    return list(array)

def compute_statistics(data):
    sorted_data = sorted(data)
    n = len(sorted_data)
    
    # Compute median
    if n % 2 == 0:
        median = (sorted_data[n // 2 - 1] + sorted_data[n // 2]) / 2
    else:
        median = sorted_data[n // 2]
    
    # Compute lower quartile
    if n % 4 == 0:
        lower_quartile = (sorted_data[n // 4 - 1] + sorted_data[n // 4]) / 2
    else:
        lower_quartile = sorted_data[n // 4]
    
    # Compute upper quartile
    if n % 4 == 0:
        upper_quartile = (sorted_data[3 * n // 4 - 1] + sorted_data[3 * n // 4]) / 2
    else:
        upper_quartile = sorted_data[3 * n // 4]

    # Compute average and variance
    np_sorted_data = np.array(sorted_data)
    average = np.mean(np_sorted_data)
    variance = np.var(np_sorted_data)

    index = int(n * 0.95)
    percentile_95 = sorted_data[index]
    
    max=sorted_data[-1]
    return {"median": median, "lower_quartile": lower_quartile, "upper_quartile": upper_quartile, "percentile_95": percentile_95, "max": max, "avg":average, "var":variance}


def profiles_multiple_stats(profiles_multiple):
    # Create an empty dictionary to store the cumulative sum of profile heights
    profile_matrix = {}

    # Iterate over each profile in profiles_multiple
    for profile in profiles_multiple:
        for stability_cutoff, heights in profile.items():
            # Iterate over each genomic coordinate and height in the profile
            for coordinate, height in heights.items():
                profile_matrix.setdefault((stability_cutoff, coordinate), [])
                profile_matrix[(stability_cutoff, coordinate)].append(height)

    profile_statistics={}

    #key = (stability_cutoff, coordinate)
    for key, heights in profile_matrix.items():
        profile_statistics[key] = compute_statistics(heights)
    
    profile_statistics_reformat={}
    for key, statistics in profile_statistics.items():
        (stability_cutoff, coordinate) = key
        profile_statistics_reformat.setdefault(stability_cutoff, {})
        profile_statistics_reformat[stability_cutoff][coordinate]=statistics
    
    return profile_statistics_reformat


def compress_random_profile_single_bin(data, statistics):
    max_position = max([int(k) for k in data.keys()])
    compressed = []
    last_stat = None
    for stat in statistics: #for each statistics
        for_stat = []
        last_value = None; last_len = 0
        for key in range(max_position+1): #for each position
            #key = str(key)
            if (key in data and data[key][stat]):
                value = data[key][stat]
            else:
                value = None
            if (last_value != value):
                if (last_len > 1):
                    for_stat.append([last_value, last_len])
                elif (last_len == 1):
                    for_stat.append([last_value])
                last_value = value; last_len = 1
            else:
                last_len += 1 

        if (last_len > 1):
            for_stat.append([last_value, last_len])
        elif (last_len == 1):
            for_stat.append([last_value])
        if (for_stat != last_stat):
            compressed.append(for_stat)
            last_stat = for_stat
        else:
            compressed.append([])
    return compressed

def compress_profile_random(data):
    statistics = ["median","lower_quartile","upper_quartile","percentile_95", "max", "avg", "var"]
    compressed = dict()
    compressed["statistics"] = statistics
    compressed["data"] = dict()
    for key in data.keys():
        compressed["data"][key] = compress_random_profile_single_bin(data[key], statistics)
    return compressed

def main():
    parser = argparse.ArgumentParser(description="read tpx stability file and compute tfo profiles")
    parser.add_argument('tpx_files', nargs='*', metavar='file', help='Input file(s)')
    parser.add_argument('-m', '--multiple_input', action='store_true', help='Multiple intput files, as in randomizations')
    parser.add_argument('-j', '--json', action='store_true', help='output json and not msgpack')
    args = parser.parse_args()

    # Check if stdin flag is provided or no files are provided
    if not args.multiple_input:
        data, length = get_TFO_profile_allSparse(sys.stdin)

        #data = json.load(sys.stdin)
        profiles = profile2ranges(data["profiles"])
        best_stability = best_stability_to_array(data["best_stability"], length)
        to_export = {"profiles": profiles, "best_stability": best_stability}

    else:
        profiles_multiple=[]
        for file_name in args.tpx_files:
            with open(file_name) as my_file:
                data, length = get_TFO_profile_allSparse(my_file)
            profiles_multiple.append(dict(data["profiles"]))
        stats = profiles_multiple_stats(profiles_multiple)
        """with open("temp_.json", "w") as o:
            o.write(json.dumps(stats))"""
        #Add compression to random profile
        """
        Uncompr: thr: { position: {"median": 1.0, "lower_quartile": 1.0, "upper_quartile": 1.0, "percentile_95": 1, "max": 1} }
        Compr: thr: [ [value, (len)] for each statistics ]
            >If (len)==1, implicit: [value]
            >For each statistic, if same as previous, empty array
                [[median...], [], [], []] if all same values
        """
        to_export = compress_profile_random(stats)

    if args.json:
        sys.stdout.write(json.dumps(to_export))
    else:
        packed_matrix = msgpack.packb(to_export, use_bin_type=True)
        sys.stdout.buffer.write(packed_matrix)        

if __name__=="__main__":
    main()
    
