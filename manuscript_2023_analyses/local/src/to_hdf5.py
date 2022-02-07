#!/usr/bin/env python3

from sys import stdin, stderr
from argparse import ArgumentParser
#from optparse import OptionParser
import pandas as pd


def main():
	usage = '''
		%prog [--usecols, --col_names, --complevel] < path/to/raw.tpx.stability.gz  > path/to/raw.tpx.stability.hdf5.gz
	\n
	.META: path/to/raw.tpx.stability.gz
	1       ssRNA
	2       TFO_start
	3       TFO_end
	4       Duplex_ID
	5       TTS_start
	6       TTS_end
	7       Score
	8       Error_rate
	9       Errors
	10      Motif
	11      Strand
	12      Orientation
	13      Guanine_rate
	14      Stability	
	
	'''

	parser = ArgumentParser(description = usage)
	
	parser.add_argument('-u', '--usecols', dest='usecols', type=int, nargs='+', default=None, help='Return a subset of the columns; if list-like, all elements must either be positional (i.e. integer indices into the document columns) or strings that correspond to column names [default: None]')
	parser.add_argument('-n', '--names', dest='names', nargs='+', default=None, help='List of column names to use [default: None]')
	parser.add_argument('-c', '--complevel', dest='complevel', type=int, default=8, help='Specifies a compression level for data. A value of 0 or None disables compression [default: 8]')
	parser.add_argument('--header', dest='header', action='store_true', default=None, help='If specified the first line will be considered the table header [default: None]')

	args = parser.parse_args()
	
	#if len(args) != 0:
	#	exit('Unexpected argument number.')
	
	# read arguments
	usecols = args.usecols	
	names = args.names
	complevel = args.complevel # add check for 0 < complevel < 9
	header = args.header

	# read data from STDIN and store in pandas dataframe
	if(names != None):
		df = pd.read_table(stdin, usecols = usecols, index_col = False, header = header, names = names)
	else:
		df = pd.read_table(stdin, usecols = usecols, index_col = False, header = header)
	
	# write the contained data to an HDF5 file using HDFStore
	df.to_hdf('/dev/stdout', key = 'group_1', complevel = complevel, mode = 'w')


if __name__ == '__main__':
	main()

