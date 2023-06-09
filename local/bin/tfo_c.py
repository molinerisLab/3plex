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
                
    return {"profiles": profiles, "best_stability": best_stability}

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
        print(">%s" % stability)
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

def main():
    get_TFO_profile_allSparse()

    path = sys.argv[1]
    #data = json.load(sys.stdin)
    profiles = profile2ranges(data["profiles"])
    best_stability = best_stability_to_range(data["best_stability"])
    to_export = {"profiles": profiles, "best_stability": best_stability}
    output_dir = os.path.join(path, "profile_range.msgpack")
    packed_matrix = msgpack.packb(to_export, use_bin_type=True)
    with open(output_dir, 'wb') as file:
        file.write(packed_matrix)

if __name__=="__main__":
    main()
    