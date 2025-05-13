#!/bin/bash

tpx_stability=$1
start=$2
end=$3
# Define the threshold value
threshold=$4

tabix $tpx_stability ssRNA:$start-$end | awk -v threshold="$threshold" -F'\t' '$14 > threshold' 