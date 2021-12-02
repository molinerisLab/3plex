#!/usr/bin/env python3
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser
#from random import shuffle
#from subprocess import Popen, PIPE
#from collections import defaultdict
#from vfork.io.colreader import Reader

A=0
C=1
G=2
T=3
N=4

PARALLEL=1
ANTIPARALLEL=0

stab_table=[[None,None,None,None,None],[None,None,None,None,None]]
#                                 A     C     G     T    N
stab_table[PARALLEL][A]=	[0.0,  0.0,  2.8,  0.0,  0.0]
stab_table[PARALLEL][C]=	[0.0,  4.5,  2.2,  2.4,  0.0]
stab_table[PARALLEL][G]=	[0.0,  2.4,  0.0,  2.6,  0.0]
stab_table[PARALLEL][T]=	[0.0,  0.0,  0.0,  3.7,  0.0]
stab_table[PARALLEL][N]=	[0.0,  0.0,  0.0,  0.0,  0.0]
stab_table[ANTIPARALLEL][A]=	[0.0,  0.0,  1.0,  0.0,  0.0]
stab_table[ANTIPARALLEL][C]=	[1.0,  3.0,  3.0,  0.0,  0.0]
stab_table[ANTIPARALLEL][G]=	[0.0,  1.0,  0.0,  2.0,  0.0]
stab_table[ANTIPARALLEL][T]=	[3.0,  1.0,  0.0,  3.5,  0.0]
stab_table[ANTIPARALLEL][N]=	[0.0,  0.0,  0.0,  0.0,  0.0]



def main():
	usage = '''
		%prog < tpx_aln > tpx_with_stab

.META: tpx_aln
	1	equence_ID
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	aln1
	15	aln2
	16	aln3
	17	aln4

.META: tpx_aln
	1	equence_ID
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	stability
	15	TFO
	16	TTS

	'''
	parser = OptionParser(usage=usage)
	
	parser.add_option('-c', '--cutoff', type=float, dest='cutoff', default=0.05, help='some help CUTOFF [default: %default]', metavar='CUTOFF')
	parser.add_option('-p', '--dump-params', dest='params_file', help='some help FILE [default: %default]', metavar='FILE')
	parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='some help [default: %default]')

	options, args = parser.parse_args()
	
	if len(args) != 0:
		exit('Unexpected argument number.')
	
	d={"A":0, "C":1, "G":2, "T":3, "N":4}
	
	for line in stdin:
		tokens = line.rstrip().split('\t')
		(Sequence_ID,TFO_start,TFO_end,Duplex_ID,TTS_start,TTS_end,Score,Error_rate,Errors,Motif,Strand,Orientation,Guanine_rate,aln1,aln2,aln3,aln4)=tokens
		if(Strand=="-"):	
			TFO_seq=aln4[9:-4].upper()
			TTS_seq=aln1[9:-4].upper()
		else:
			TFO_seq=aln1[9:-4].upper()
			TTS_seq=aln4[9:-4].upper()
			
		P=PARALLEL
		if(Orientation=="A"):
			P=ANTIPARALLEL
		
		#print(TFO_seq, TTS_seq)
		stability=0
		for (i,c) in enumerate(TTS_seq):
			#print(c,TFO_seq[i])
			#stability+=stab_table[P][c][TFO_seq[i]]
			stability+=stab_table[P][d[c]][d[TFO_seq[i]]]
		print("\t".join(tokens[:-4]+[str(round(stability,1)),aln1,aln2,aln3,aln4]))
		#print("\t".join(tokens[:-4]+[str(round(stability,1)),TFO_seq,TTS_seq]))	

if __name__ == '__main__':
	main()

