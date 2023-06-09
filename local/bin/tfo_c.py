#!/bin/env python
#import pandas as pd
#import numpy as np
#from scipy.sparse import dok_matrix
from collections import defaultdict
import json
import sys


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
                
    print(json.dumps({"profiles": profiles, "best_stability": best_stability}))


if __name__=="__main__":
    get_TFO_profile_allSparse()