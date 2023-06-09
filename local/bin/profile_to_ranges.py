#!/bin/env python
#import pandas as pd
#import numpy as np
#from scipy.sparse import dok_matrix
import msgpack
import json
import sys
import os 

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


if __name__ == '__main__':
    path = sys.argv[1]
    data = json.load(sys.stdin)
    profiles = profile2ranges(data["profiles"])
    best_stability = best_stability_to_range(data["best_stability"])
    to_export = {"profiles": profiles, "best_stability": best_stability}

    output_dir = os.path.join(path, "profile_range.msgpack")
    """with open('stability_big.cut_sort.profile_range.json', 'w') as file:
        file.write(json.dumps(to_export))"""
    
    packed_matrix = msgpack.packb(to_export, use_bin_type=True)
    with open(output_dir, 'wb') as file:
        file.write(packed_matrix)

