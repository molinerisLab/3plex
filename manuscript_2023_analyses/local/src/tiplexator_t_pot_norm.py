#!/usr/bin/env python3.6
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser
#from random import shuffle
#from subprocess import Popen, PIPE
#from collections import defaultdict
#from vfork.io.colreader import Reader



def main():
	usage = '''
		%prog < STDIN
.META: stdin
	1	ssRNA	TERC
	2	duplex_ID chirp_peak_1000;TERC
	3	tpx_count	9
	4	duplexlength	652
	5	oligolength	541
	'''

	parser = OptionParser(usage=usage)
	
	parser.add_option('-l', '--minLength', type=int, dest='minLength', default=10, help='[default: %default]', metavar='CUTOFF')
	parser.add_option('-L', '--maxLength', type=int, dest='maxLength', default=30, help='-1 to no limit [default: %default]', metavar='CUTOFF')

	options, args = parser.parse_args()
	
	if len(args) != 0:
		exit('Unexpected argument number.')
	
	for line in stdin:
		(ssRNA,duplex_ID,tpx_count,duplexlength,oligolength) = line.rstrip().split('\t')
		#oligolength===ssRNA_length
		duplexlength=int(duplexlength)
		oligolength=int(oligolength)
		tpx_count=int(tpx_count)
		
		if options.maxLength==-1:
			maxLength = min((oligolength,duplexlength))
		else:
			maxLength = min((options.maxLength,oligolength,duplexlength))

		n=0
		for i in range(options.minLength,maxLength+1):
  			n += (duplexlength-i+1)*(oligolength-i+1)
		n*=2;
	
		print("\t".join((str(i) for i in (ssRNA,duplex_ID,tpx_count,duplexlength,oligolength,n,float(tpx_count)/n))))

if __name__ == '__main__':
	main()

