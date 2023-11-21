#!/usr/bin/env python3

from __future__ import with_statement
import sys
import numpy as np
import argparse


def main():

	parser = argparse.ArgumentParser(
		description="Compute the modified z-score of each RNA sequence position from a lunp input file (RNAplfold output).",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)

	args = parser.parse_args()
	
	# read stdin and put values in array
	score_lst = [ float(val.strip()) for val in sys.stdin if val != "NA\n" ]
	score_nparr = np.array(score_lst)

	# compute median and median absolute deviation from the median 
	median = np.median(score_nparr)
	mad = np.median(np.absolute(score_nparr - median))

	# modified z scores array
	modif_zscore = (score_nparr - median) / mad

        # print stdout
	for i in modif_zscore:
		print("{:.4f}".format(i))

if __name__ == '__main__':
        main()
