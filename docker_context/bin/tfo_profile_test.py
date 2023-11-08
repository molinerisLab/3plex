#!/bin/env python
#import pandas as pd
#import numpy as np
#from scipy.sparse import dok_matrix
from collections import defaultdict
import json
import sys
import msgpack
import sys
import os 
import numpy as np

def get_TFO_profile_allSparse():
    # stdin has 3 columns:
    #
    #  1       TFO_start
    #  2       TFO_end
    #  3       Stability
    #
    # without header
    #
    # assume data sorted on desending Stability

    #profiles=dok_matrix((stability_max, ssRNA_len),dtype=np.uintc)
    profiles = {}
    profile_current = defaultdict(lambda: 0)
    stability_pre=None
    best_stability = {}
    max_len = 0
    
    for line in sys.stdin:
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

def ranges_back_to_dict(ranges):
    #profiles: dict pos - n
    #ranges: [ n, pos ] or [n, pos_middle, len]
    #profiles: for stability in profiles.keys(): profile = profiles[stability]
    #ranges: same thing
    prof = dict()
    for stability in ranges.keys():
        prof[stability] = defaultdict(lambda: 0)
        profiles = ranges[stability]
        for profile in profiles:
            if (len(profile)==2):
                prof[stability][profile[1]] = profile[0]
            elif(len(profile)==3):
                min_ = profile[1] - ((profile[2]-1)/2); max_ = profile[1] + ((profile[2]-1)/2);
                min_ = int(min_); max_ = int(max_)
                #print(profile)
                for i in range(min_, max_+1):
                    prof[stability][i] = profile[0]
                    #print(f"{i} - ")
                #exit()
            else:
                print(profile); exit()
    return prof


def main():
    data, length = get_TFO_profile_allSparse()
    profile_ranges = profile2ranges(data["profiles"])
    data = data["profiles"]
    data_r = ranges_back_to_dict(profile_ranges)
    #print(data_r); exit()
    n_v = 0; n_err = 0;
    for stability in data.keys():
        profile = data[stability]
        for value in profile.keys():
            n_v += 1
            if (profile[value] != data_r[stability][value]):
                n_err += 1
                print(f"{value}: - {profile[value]} != {data_r[stability][value]}")
    print(n_v); print(n_err)

if __name__=="__main__":
    main()
    