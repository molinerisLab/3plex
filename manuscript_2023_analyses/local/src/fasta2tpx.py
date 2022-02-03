#!/usr/bin/env python
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser
#from random import shuffle
#from subprocess import Popen, PIPE
#from collections import defaultdict
from Bio import SearchIO


def main():
	usage = '''
		%prog PARAM1 PARAM2 < STDIN
	'''
	parser = OptionParser(usage=usage)
	
	parser.add_option('-c', '--cutoff', type=float, dest='cutoff', default=0.05, help='some help CUTOFF [default: %default]', metavar='CUTOFF')
	parser.add_option('-p', '--dump-params', dest='params_file', help='some help FILE [default: %default]', metavar='FILE')
	parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='some help [default: %default]')

	options, args = parser.parse_args()
	
	if len(args) != 0:
		exit('Unexpected argument number.')

	hits = SearchIO.read(stdin,"fasta-m10")
	
	for hit in hits:
		for hps in hit.hsps:
#			print(hps)
#			print(">{}\t{}\t{}\t{}\t{}".format(hps.query_id, hps.hit_id, hps.bitscore, hps.evalue, hps.z_score))
			for fragment in hps.fragments:
				query_id	= fragment.query_id
				target_id	= fragment.hit_id
				query_start	= fragment.query_start
				query_end	= fragment.query_end
				target_start	= fragment.hit_start
				target_end	= fragment.hit_end
				query_strand	= fragment.query_strand
				target_strand	= fragment.hit_strand
				query_fragment	= fragment.query.seq
				target_fragment = fragment.hit.seq
				similarity 	= fragment.aln_annotation['similarity']
			
				print("\t".join((str(s) for s in (
					query_id,
					target_id,
					query_start,
					query_end,	
					target_start,
					target_end,
					query_strand,
					target_strand,	
					query_fragment,	
					target_fragment, 
					similarity
				))))


	
if __name__ == '__main__':
	main()

