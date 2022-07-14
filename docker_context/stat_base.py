#!/usr/bin/env python
#
# Copyright 2008,2010 Gabriele Sales <gbrsales@gmail.com>, 2009,2010 Ivan Molineris <ivan.molineris@gmail.com>

from __future__ import division
from itertools import groupby
from operator import itemgetter 
from optparse import OptionParser
from sys import stdin
from math import sqrt
from vfork.io.colreader import Reader
from vfork.io.util import safe_rstrip
from vfork.util import exit, format_usage
from operator import mul
from math import ceil, floor
from numpy import std,mean


class ValueReader(object):
	def __init__(self, fd):
		self.fd = fd
		self.lineno = 0
	
	def iter_groups(self, allready_sorted):
		values = []
		for line in self.fd:
			self.lineno += 1
			values.append(self._parse_value(safe_rstrip(line)))
		return [(None, values)]
	
	def _parse_value(self, repr):
		try:
			return float(repr)
		except ValueError:
			exit('Invalid value at line %d (%s).' % (self.lineno, repr))

class GroupReader(object):
	def __init__(self, fd, sorted):
		self.fd = fd
		self.sorted = sorted

	def iter_groups(self, allready_sorted):
		last_grp = None
		for grp, lines in groupby(self._iter_lines(), itemgetter(0)):
			if not self.sorted and last_grp is not None and grp < last_grp:
				exit('Disorder found in input while processing group: ' + grp)
			last_grp = grp

			values = [ l[1] for l in lines ]
			yield grp, values

	def _iter_lines(self):
		reader = Reader(self.fd, '0s,1f', False)
		for k, v in reader:
			yield k, v

def stdev(values, mean):
	if len(values) > 2:
		c = 0
		for v in values:
			c += (v-mean)**2 # -2**2 = -4 ma solo per questione di precedenza degli operatori invece a=-2; a**2; a==4
		stdev = sqrt(c / (len(values)-1))
	else:
		stdev = None
	return stdev

def skewness(histogram_binned_values):# pearson_momet_coefficient taken from geneBody_coverage.py in rseqc
	'''measure skewness'''
	mid_value = histogram_binned_values[int(len(histogram_binned_values)/2)]
	sigma = std(histogram_binned_values, ddof=1)
	tmp = []
	for i in histogram_binned_values:
		tmp.append(((i - mid_value)/sigma)**3)
	return mean(tmp)

def label(val,label):
	global options
	if val == None:
		val = options.null_value
	if options.labels: 
		return label + ":" + str(val)
	else: 
		return str(val)

def main():
	global options
	parser = OptionParser(usage=format_usage('''
		%prog [OPTIONS] <VALUES

		Computes basic statistic parameters for the input values.

		.META stdout
			1	grp
			2	mean
			3	median
			4	stdev_mean
			5	min
			6	max
			7	len
			8	sum
			9	prod

		NOTE: the order ot the output columns is preserved but only the column associated with given flags are printed
	'''))
	parser.add_option('-a', '--mean', dest='mean', action='store_true', default=False, help='compute the mean value')
	parser.add_option('-m', '--median', dest='median', action='store_true', default=False, help='compute the median value')
	parser.add_option('-q', '--percentile', dest='percentile', type='float', default=None, help='compute the percentile Q (in fraction, i.e. the median is -q 0.5)', metavar="Q")
	parser.add_option('-Q', '--percentile2', dest='percentile2', type='float', default=None, help='compute the percentile Q2 (when you want to compute twho different percentile in the same run)', metavar="Q")
	parser.add_option('-s', '--stdev', dest='stdev', action='store_true', default=False, help='compute the standard deviation')
	parser.add_option('-S', '--mean_stdev', dest='mean_stdev', action='store_true', default=False, help='compute the standard deviation of the mean (stdev/sqrt(N-1))')
	parser.add_option('-l', '--min', dest='min', action='store_true', default=False, help='compute the min')
	parser.add_option('-b', '--max', dest='max', action='store_true', default=False, help='compute the max')
	parser.add_option('-c', '--count', dest='count', action='store_true', default=False, help='count the number of rows')
	parser.add_option('-t', '--tot', dest='summ', action='store_true', default=False, help='summ all the values')
	parser.add_option('-p', '--prod', dest='prod', action='store_true', default=False, help='multiply all the values')
	parser.add_option('-k', '--skewness', dest='skewness', action='store_true', default=False, help='compute skewness (pearson moment coefficient) assuming as input ordered hisogram binned values (like geneBody_coverage)')

	parser.add_option('-z', '--nonzero', dest='nonzero', action='store_true', default=False, help='count occurrences of non zero values')
	parser.add_option('-g', '--groups', dest='groups', action='store_true', default=False, help='group values by the key in the first column')
	parser.add_option('-o', '--sorted', dest='sorted', action='store_true', default=False, help='assume the input is sorted on the first column')
	parser.add_option('-n', '--ignore-empty', dest='ignore_empty', action='store_true', default=False, help='do not exit with an error when there are no records')
	parser.add_option('-N', '--null_value', dest='null_value', type='string', default="", help='use ND insead of "" when a value is not defined (like the stdev of only one data) ', metavar="ND")
	parser.add_option('-L', '--labels', dest='labels', action='store_true', default=False, help='prepend each value reported in stdout with a label describing what it is) ')
	options, args = parser.parse_args()
	
	if len(args) != 0:
		exit('Unexpected argument number.')
	
	if options.groups:
		reader = GroupReader(stdin, options.sorted)
	else:
		reader = ValueReader(stdin)
	
	s = None
	for grp, values in reader.iter_groups(options.sorted):
		out = []
		s = None
		
		if options.groups:
			out.append(grp)
		
		if options.mean or options.stdev:
			mean = sum(values) / len(values)
			if options.mean:
				out.append(label(mean,"mean"))
		
		if options.median:
			values.sort()
			if len(values) % 2 == 0:
				idx = len(values) // 2
				median = (values[idx] + values[idx-1]) / 2
			else:
				median = values[len(values)//2]
			out.append(label(median,"median"))

		if options.percentile is not None:
			if not options.median:
				values.sort()
			percentile = len(values)*options.percentile
			if percentile % 1 != 0:
				percentile = int(ceil(percentile)) #For a finite population of N values indexed 1,...,N from lowest to highest, the kth q-quantile of this population can be computed via the value of I_p = N \frac{k}{q}. If I_p is not an integer, then round up to the next integer to get the appropriate index;http://en.wikipedia.org/wiki/Quantile#Quantiles_of_a_population
				percentile-=1 #For a finite population of N values indexed 1,...,N not from 0
				percentile_val = values[percentile]
			else:
				percentile = int(percentile) -1
				percentile_val = (values[percentile] + values[percentile + 1]) / 2 # if I_p is an integer then any number from the data value at that index to the data value of the next can be taken as the quantile, and it is conventional (though arbitrary) to take the average of those two values (see 
			out.append(label(percentile_val,"percentile"))
		if options.percentile2 is not None:
			if not options.median and not options.percentile:
				values.sort()
			percentile = len(values)*options.percentile2
			if percentile % 1 != 0:
				percentile = int(ceil(percentile)) #For a finite population of N values indexed 1,...,N from lowest to highest, the kth q-quantile of this population can be computed via the value of I_p = N \frac{k}{q}. If I_p is not an integer, then round up to the next integer to get the appropriate index;http://en.wikipedia.org/wiki/Quantile#Quantiles_of_a_population
				percentile-=1 #For a finite population of N values indexed 1,...,N not from 0
				percentile_val = values[percentile]
			else:
				percentile = int(percentile) -1
				percentile_val = (values[percentile] + values[percentile + 1]) / 2 # if I_p is an integer then any number from the data value at that index to the data value of the next can be taken as the quantile, and it is conventional (though arbitrary) to take the average of those two values (see 
			out.append(label(percentile_val,"percentile2"))
		
		if options.stdev:
			s = stdev(values, mean)
			out.append(label(s,"stdev"))

		if options.mean_stdev:
			if len(values) > 1:
				if s is None:
					s = stdev(values, mean)
				if s is None:
					out.append(label(None,"mean_stdev"))
				else:
					out.append(label( s / sqrt(len(values) - 1),"mean_stdev") )
			else:
				out.append(label(None,"mean_stdev"))

		if options.min:
			out.append(label(min(values),"min"))
		
		if options.max:
			out.append(label(max(values),"max"))
		
		if options.count:
			out.append(label(len(values),"len"))

		if options.summ:
			out.append(label(sum(values),"sum"))

		if options.nonzero:
			out.append(label(sum(1 for v in values if v!=0.),"nonzero"))

		if options.prod:
			out.append(label(reduce(mul, values),"prod"))

                if options.skewness:
                        out.append(label(skewness(values),"skewness"))

		print '\t'.join(out)

if __name__ == '__main__':
	main()

