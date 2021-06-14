#!/usr/bin/env python3
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser
#from random import shuffle
#from subprocess import Popen, PIPE
from collections import defaultdict
from itertools import product
#from vfork.io.colreader import Reader

def main():
	usage = '''
%prog < STDIN
'''
	parser = OptionParser(usage=usage)
	
	#parser.add_option('-c', '--cutoff', type=float, dest='cutoff', default=0.05, help='some help CUTOFF [default: %default]', metavar='CUTOFF')
	#parser.add_option('-p', '--dump-params', dest='params_file', help='some help FILE [default: %default]', metavar='FILE')
	#parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='some help [default: %default]')

	options, args = parser.parse_args()
	
	if len(args) != 0:
		exit('Unexpected argument number.')
	
	data = defaultdict(lambda: defaultdict(float))
	d_sets=defaultdict(set)
	scores=set()

	for line in stdin:
		k1,k2,v = line.rstrip().split('\t')
		v=float(v)
		data[k1][k2]=v
		scores.add(v)

	for s in sorted(scores):
		d_sets=defaultdict(set)
		for k1 in data.keys():
			for k2 in data[k1].keys():
				if data[k1][k2]>=s:
					d_sets[k1].add(k2)
		
		for k1_1,k1_2 in product(d_sets.keys(),repeat = 2):
			u = len(d_sets[k1_1] | d_sets[k1_2])
			i = len(d_sets[k1_1] & d_sets[k1_2])
			jaccard = i/u
			print("%s\t%s\t%f\t%d\t%u\t%f" % (k1_1,k1_2,s,i,u,jaccard))

if __name__ == '__main__':
	main()

