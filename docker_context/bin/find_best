#!/usr/bin/env python
#
# Copyright 2011,2014 Gabriele Sales <gbrsales@gmail.com>

from collections import defaultdict
from itertools import groupby
from operator import itemgetter
from optparse import OptionParser
from sys import stdin
#from vfork.io.util import safe_rstrip
#from vfork.util import format_usage, exit, ignore_broken_pipe

SKIP_HEADER=False


def main():
    parser = OptionParser(usage='''
      Usage: %prog [OPTIONS] ID_COL[:ID_COL...] SCORE_COL[:SCORE_COL...] <TABLE

      For each group of rows with the same value in ID_COL(s), select the row
      with the highest SCORE.
    ''')
    parser.add_option('-s', '--sorted', dest='sorted', default=False, action='store_true', help='assume the input is already sorted by ID')
    parser.add_option('-r', '--reverse', dest='reverse', default=False, action='store_true', help='take the lowest SCORE(s)')
    parser.add_option('-m', '--multiple', dest='multiple', default=False, action='store_true', help='in case of ties, select all the rows with higher SCORE')
    parser.add_option('-n', '--many', dest='many', type='int', help='select N rows', metavar='N')
    parser.add_option('-a', '--absolute', dest='absolute', default=False, action='store_true', help='sort by absolute values.')
    parser.add_option('-H', '--skip_header', dest='skip_header', default=False, action='store_true', help='skip the first row and report it as is in stdout')

    options, args = parser.parse_args()

    if len(args) != 2:
        exit('Unexpected argument number.')
    elif options.many is not None:
        if options.multiple:
            exit('Invalid options: "--multiple" cannot be used with "--many".')
        elif options.many < 1:
            exit('Invalid value of N: %d' % options.many)

    global SKIP_HEADER
    SKIP_HEADER=options.skip_header

    it = iter_sorted_records if options.sorted else iter_unsorted_records
    id_cols = ColumnGroup(args[0], 'ID_COL(s)')
    score_cols = ScoreColumns(args[1], 'SCORE_COL(s)', options.absolute)
    get_first = itemgetter(0)

    for id, grp in it(stdin, id_cols, score_cols):
        grp.sort(key=get_first, reverse=not options.reverse)

        if options.many is not None:
            grp = grp[:options.many]

        elif options.multiple:
            best_score = grp[0][0]
            top = 0
            for i in range(1, len(grp)):
                if grp[i][0] == best_score:
                    top = i
                else:
                    break
            grp = grp[:top+1]

        else:
            grp = grp[:1]

        for _, line in grp:
            print(line)


class ColumnGroup(object):

    def __init__(self, str, descr):
        try:
            self.idxs = [ int(i)-1 for i in str.split(':') ]
            for idx in self.idxs:
                if idx < 0: raise ValueError
        except ValueError:
            exit('Invalid %s: %s' % (descr, str))

    def select(self, tokens):
        return tuple(tokens[i] for i in self.idxs)

class ScoreColumns(ColumnGroup):

    def __init__(self, str, descr, absolute):
        super(ScoreColumns, self).__init__(str, descr)

        if absolute:
            self.conv = lambda x: abs(float(x))
        else:
            self.conv = float

    def select(self, tokens):
        vs = super(ScoreColumns, self).select(tokens)
        conv = self.conv
        return tuple(conv(v) for v in vs)


def iter_records(fd, id_cols, score_cols):
    for lineno, line in enumerate(stdin):
        line=line.rstrip("\r\n")
        if SKIP_HEADER and lineno==0:
            print(line)
            continue
        
        tokens = line.split('\t')

        try:
            id = id_cols.select(tokens)
            score = score_cols.select(tokens)

        except IndexError:
            exit('Insufficient columns at line %d: %s' % (lineno+1, line))
        except ValueError:
            exit('Invalid score value (must be numeric) at line %d: %s' % (lineno+1, line))

        yield id, (score, line)

def iter_sorted_records(fd, id_cols, score_cols):
    for id, grp in groupby(iter_records(fd, id_cols, score_cols), itemgetter(0)):
        yield id, [ e[1] for e in grp ]

def iter_unsorted_records(fd, id_cols, score_cols):
    store = defaultdict(list)
    for id, info in iter_records(fd, id_cols, score_cols):
        store[id].append(info)

    for id in sorted(store):
        yield id, store[id]


if __name__ == '__main__':
    main()
