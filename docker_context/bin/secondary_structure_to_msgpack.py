#!/bin/env python
import sys
import msgpack
import os 

def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[:-len(suffix)]
    return input_string

def main():
    output_path = sys.argv[1]
    signal = []
    for line in sys.stdin:
        line = remove_suffix(line,"\n")
        signal.append(line)  
    packed_signal = msgpack.packb(signal, use_bin_type=True)
    with open(output_path, 'wb') as file:
        file.write(packed_signal)

if __name__=="__main__":
    main()
    