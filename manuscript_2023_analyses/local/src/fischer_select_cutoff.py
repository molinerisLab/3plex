#!/usr/bin/env python3
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser
from scipy.stats import fisher_exact
#from random import shuffle
#from subprocess import Popen, PIPE
from collections import defaultdict
#from vfork.io.colreader import Reader



def main():
	usage = '''
%prog < STDIN


	.META: STDIN
		1	key	a test is calcualted for each key separately
		2	score	numeric
		3	pos_neg	eiter in {"pos","neg"} or in {0,1}

	for each possible cutoff in scores, compute a contingency table of greater / lower / positive / negativa and apply a fisher exact test

	.META: STDOUT
		1	key
		2	greater and positive
		3	greater and negative
		4	lower and positive
		5	lower and negative
		6 	oddsratio
		7	pvalue

	oddsratio is the ratio of greater/lower positive cases, normalized over the same ratio in negative

	The calculated odds ratio is different from the one R uses. This scipy implementation returns the (more common) "unconditional Maximum Likelihood Estimate", while R uses the "conditional Maximum Likelihood Estimate".
'''

	parser = OptionParser(usage=usage)
	
	parser.add_option('-a', '--alternative', type=str, dest='alternative', default="two-sided", help='Defines the alternative hypothesis, The following options are available: two-sided, less, greater [default: %default]', metavar='SIDE')
	#parser.add_option('-p', '--dump-params', dest='params_file', help='some help FILE [default: %default]', metavar='FILE')
	#parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='some help [default: %default]')

	options, args = parser.parse_args()
	
	if len(args) != 0:
		exit('Unexpected argument number.')
	
	#for id, sample, raw, norm in Reader(stdin, '0u,1s,2i,3f', False):
	data=defaultdict(lambda: defaultdict(list))
	#scores_min=dict()
	#scores_max=dict()
	for line in stdin:
		key,score,pos_neg = line.rstrip().split('\t')
		score=float(score)
		if pos_neg=="pos":
			pos_neg=1
		if pos_neg=="neg":
			pos_neg=0
		pos_neg=int(pos_neg)
		assert(pos_neg==1 or pos_neg==0)

		data[key][pos_neg].append(score)

		#s=scores_min.get(key)
		#if s==None or s>score:
		#	scores_min[key]=s
		#s=scores_max.get(key)
		#if s==None or s<score:
		#	scores_max[key]=s
	
	
	for k in data.keys():
		scores = data[k][0]+data[k][1]
		for s in sorted(set(scores)):
			g_p = len([i for i in data[k][1] if i>=s])
			g_n = len([i for i in data[k][0] if i>=s])
			l_p = len([i for i in data[k][1] if i<s])
			l_n = len([i for i in data[k][0] if i<s])
			oddsratio, pvalue = fisher_exact([[g_p,g_n],[l_p,l_n]], options.alternative)

			print("%s\t%s\t%d\t%d\t%d\t%d\t%f\t%g" % (k, s, g_p, l_p, g_n, l_n, oddsratio, pvalue))

if __name__ == '__main__':
	main()

