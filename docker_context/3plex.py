#!/usr/bin/env python3
from __future__ import with_statement

import os
import re
import argparse
import yaml
import tempfile
import subprocess

def execute(cmd, cwd):
    print(cmd)
    #with subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True, bufsize=1) as p:
    with subprocess.Popen(['/bin/bash', '-c', cmd], cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, bufsize=1) as p:
        for line in p.stdout:
            print(line, end='')
        for line in p.stderr:
            print(line, end='')
        return_code = p.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

def get_ssRNA_seq_name(ssRNA):
	with open(ssRNA) as fasta:
		header=fasta.readline()
		m=re.match(r'>([^\s]+)',header)
		return m.group(1)


parser = argparse.ArgumentParser(
    description="Given an RNA sequence and a list of DNA sequences, compute all possible triplexes that satisfy the constraints associating a thermal stability score. RNA secondary structure prediction can be used to exclude RNA nucleotides from search.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    epilog="Have a good triplex journey!")

parser.add_argument("ssRNA", metavar="ssRNA.fa", help="The RNA sequence in fasta format. The file must contain only one sequence.")
parser.add_argument("dsDNA", metavar="dsDNA.fa", help="The DNA seeuqnces. Multiple fasta formatted sequences, intende to be double stranded DNA potentially bound by the RNA sequence.")
parser.add_argument("out_dir", metavar="out_dir", help="The output directory. Out_dir should already exist. Pay attention to mount appropraite volumes if you ar using docker (e.g. `docker run --rm -v $(PWD):$(PWD) 3plex:v0.0.1 ssRNA.fa dsDNA.fa $(PWD)` sould work as expected, `docker run --rm -v $(PWD):$(PWD) 3plex:v0.0.1 ssRNA.fa dsDNA.fa .` should not)")
parser.add_argument("--bed", metavar="dsDNA.bed", dest="dsDNA_bed", type=str, default=None, help="Genomic coordiantes of the DNA sequences in bed format, the 4th column must contain the same identifiers of sequences in dsDNA.fa")
parser.add_argument("-j", "--jobs", metavar="CPUS", dest="jobs", type=int, default=1, help="Number of parallel threads.")
parser.add_argument("-l", "--min_length", metavar="N", dest="triplexator_min_length", type=int, default=10, help="Minimum triplex length required. Allowed values: N>5.")
parser.add_argument("-e", "--error_rate", metavar="E", dest="triplexator_error_rate", type=int, default=20, choices=range(100), help="Maximal percentage of error allowed in a triplex.")
parser.add_argument("-s", "--single_strandedness_cutoff", metavar="S", dest="RNAplfold_single_strandedness_cutoff", type=int, default=0, choices=range(100), help="Percentage of masked RNA nucleotides based on RNAplfold base pairing probabilities.")
parser.add_argument("-c", "--consecutive_errors", metavar="C", dest="triplexator_consecutive_errors", type=int, default=3, help="Maximum number of consecutive errors allowed in a triplex.")
parser.add_argument("-g", "--guanine_rate", metavar="G", dest="triplexator_guanine_rate", type=int, default=10, choices=range(100), help="Minimal guanine percentage required in any TTS.")
parser.add_argument("-r", "--filter_repeat", metavar="R", dest="triplexator_filter_repeat", type=bool, default=False, help="If enabled, exclude repeat and low complexity regions.")
parser.add_argument("-L", "--max_length", metavar="M", dest="triplexator_max_length", type=int, default=-1, help="Maximum triplex length permitted, M=-1 imply no limits.")
parser.add_argument("-t", "--triplexator_other_parameters", metavar="T", dest="triplexator_other_parameters", type=str, default="", help="Additional triplexator parameters passed as a sting (e.g. -t '-mamg 90 -E 4'). Triplexator output format will not change.")
parser.add_argument("--ucsc_dark_gray", metavar="G", dest="TTS_bed_ucsc_dark_gray", type=int, default=843, choices=range(1000), help="TTS bed UCSC dark gray")
parser.add_argument("--dark_gray_stability", metavar="G", dest="TTS_bed_ucsc_dark_gray_stability", type=int, default=43, help="10%% of TTS in paper.")
parser.add_argument("--RNAplfold_window_size", metavar="S", dest="RNAplfold_window_size", type=int, default=200, help="RNAplfold: average pair probabilities over windows of specified size.")
parser.add_argument("--RNAplfold_span_size", metavar="S", dest="RNAplfold_span_size", type=int, default=150, help="RNAplfold: maximum separation of a base pair permitted.")
parser.add_argument("--RNAplfold_unpaired_window", metavar="S", dest="RNAplfold_unpaired_window", type=int, default=8, help="RNAplfold: mean probability that regions of specified length are unpaired.")
parser.add_argument("--snakefile", metavar="file", dest="snakefile", type=str, default="/opt/3plex/Snakefile", help=" ")
parser.add_argument("--no_env", dest="no_env", type=bool, default=False, help="Do not load the conda environment, useful when running the script outside of the docker image.")
args = parser.parse_args()

if args.triplexator_filter_repeat:
	args.triplexator_filter_repeat="on"
else:
	args.triplexator_filter_repeat="off"

config={}
config["triplexator"]={}
config["RNAplfold"]={}
config["TTS_bed"]={}

config["triplexator"]["min_length"]=args.triplexator_min_length
config["triplexator"]["max_length"]=args.triplexator_max_length
config["triplexator"]["error_rate"]=args.triplexator_error_rate
config["triplexator"]["guanine_rate"]=args.triplexator_guanine_rate
config["triplexator"]["filter_repeat"]=args.triplexator_filter_repeat
config["triplexator"]["consecutive_errors"]=args.triplexator_consecutive_errors
config["triplexator"]["other_parameters"]=args.triplexator_other_parameters

config["RNAplfold"]["window_size"]=args.RNAplfold_window_size
config["RNAplfold"]["span_size"]=args.RNAplfold_span_size
config["RNAplfold"]["unpaired_window"]=args.RNAplfold_unpaired_window
config["RNAplfold"]["unpaired_window"]=args.RNAplfold_unpaired_window
config["RNAplfold"]["single_strandedness_cutoff"]=args.RNAplfold_single_strandedness_cutoff

config["TTS_bed"]["ucsc_dark_gray"]=args.TTS_bed_ucsc_dark_gray
config["TTS_bed"]["dark_gray_stability"]=args.TTS_bed_ucsc_dark_gray_stability

config_yaml = yaml.dump(config)


if not os.path.exists(args.out_dir):
	os.mkdir(args.out_dir)

tmpdir = tempfile.mkdtemp(prefix=args.out_dir + "/" + "3plex_tmp_directory_", dir=".")

ssRNA_name=get_ssRNA_seq_name(args.ssRNA)

os.symlink(args.ssRNA, "{tmpdir}/{ssRNA_name}.fa".format(ssRNA_name=ssRNA_name, tmpdir=tmpdir))
	
os.symlink(args.dsDNA, tmpdir+"/"+os.path.basename(args.dsDNA))
os.symlink(args.snakefile, tmpdir+"/"+os.path.basename(args.snakefile))
if args.dsDNA_bed:
	os.symlink(args.dsDNA_bed, tmpdir+"/"+os.path.basename(dsDNA_bed))

with open(tmpdir+"/config.yaml","w") as config_file:	
    config_file.write(config_yaml)

bashCommand = """
source /etc/profile; 
conda activate 3plex_v0.1;
""" 
if args.no_env:
	bashCommand=""

bashCommand+="""
snakemake --snakefile {snakefile} -j {jobs} {ssRNA}_ssmasked-{dsDNA}.tpx.summary.gz {ssRNA}_ssmasked-{dsDNA}.tpx.stability.gz && \
mv {ssRNA}_ssmasked-{dsDNA}.tpx.summary.gz {ssRNA}_ssmasked-{dsDNA}.tpx.stability.gz ../
""".format(
	jobs=args.jobs, 
	snakefile=os.path.basename(args.snakefile),
	ssRNA=ssRNA_name,
	dsDNA=os.path.splitext(os.path.basename(args.dsDNA))[0]
)
execute(bashCommand, tmpdir)

