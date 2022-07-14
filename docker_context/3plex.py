#!/usr/bin/env python3
from __future__ import with_statement

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



parser = argparse.ArgumentParser(
    description="It is a long established fact that a reader will be distracted by the readable content of a page when looking at its layout. The point of using Lorem Ipsum is that it has a more-or-less normal distribution of letters, as opposed to using 'Content here, content here', making it look like readable English. Many desktop publishing packages and web page editors now use Lorem Ipsum as their default model text, and a search for 'lorem ipsum' will uncover many web sites still in their infancy. Various versions have evolved over the years, sometimes by accident, sometimes on purpose (injected humour and the like).",
     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
     epilog="Have a good triplx journey!")

parser.add_argument("ssRNA", metavar="ssRNA.fa", help="The RNA sequence in fasta format. The file must contain only one sequence and the filename ([name].fa) must be the same name of the sequence in the fasta header (>[name]).")
parser.add_argument("dsDNA", metavar="dsDNA.fa", help="The DNA seeuqnces. A fasta file containing many sequences, intende to be double stranded DNA potentially bound by the RNA sequence.")
parser.add_argument("out_dir", metavar="out_dir", help="The output directory. Out_dir should already exists. Pay attention to mount appropraite volumes if you ar using docker (e.g. `docker run --rm -v $(PWD):$(PWD) 3plex:v0.0.1 ssRNA.fa dsDNA.fa $(PWD)` sould work as expected, `docker run --rm -v $(PWD):$(PWD) 3plex:v0.0.1 ssRNA.fa dsDNA.fa .` maybe not)")
parser.add_argument("--bed", metavar="dsDNA.bed", dest="dsDNA_bed", type=str, default=None, help="Genomic coordiantes of the DNA seeuqnces in bed format, the 4th column must contain the same identifiers of sequences in dsDNA.fa")
parser.add_argument("-j", "--jobs", metavar="CPUS", dest="jobs", type=int, default=1, help="Number of parallels threads.")
parser.add_argument("-l", "--min_length", metavar="N", dest="triplexator_min_length", type=int, default=10, help="Specifies the minimum length of triplex, allowed values: N>5.")
parser.add_argument("-e", "--error_rate", metavar="E", dest="triplexator_error_rate", type=int, default=20, choices=range(100), help="Specifies the percentage of errors allowed in a TPX.")
parser.add_argument("-s", "--single_strandedness_cutoff", metavar="S", dest="RNAplfold_single_strandedness_cutoff", type=int, default=0, choices=range(100), help=" ")
parser.add_argument("-c", "--consecutive_errors", metavar="C", dest="triplexator_consecutive_errors", type=int, default=3, help="Specifies the number of consecutive_errors allowed in a TPX.")
parser.add_argument("-g", "--guanine_rate", metavar="G", dest="triplexator_guanine_rate", type=int, default=10, choices=range(100), help="Specifies the minimum guanine percentage in any TTS.")
parser.add_argument("-r", "--filter_repeat", metavar="R", dest="triplexator_filter_repeat", type=bool, default=False, help="Filter repeat and low complexity regions.")
parser.add_argument("-L", "--max_length", metavar="M", dest="triplexator_max_length", type=int, default=-1, help="Specifies the maximum length of triplex, M=-1 imply no limits.")
parser.add_argument("-t", "--triplexator_other_parameters", metavar="T", dest="triplexator_other_parameters", type=str, default="", help="Additional triplexator parameters, as sting (e.g. -t '-mamg 90 -E 4'). Do not change the triplexator output format.")
parser.add_argument("--ucsc_dark_gray", metavar="G", dest="TTS_bed_ucsc_dark_gray", type=int, default=843, choices=range(1000), help=" ")
parser.add_argument("--dark_gray_stability", metavar="G", dest="TTS_bed_ucsc_dark_gray_stability", type=int, default=43, help="10%% of TTS in paper.")
parser.add_argument("--RNAplfold_window_size", metavar="S", dest="RNAplfold_window_size", type=int, default=200, help=" ")
parser.add_argument("--RNAplfold_span_size", metavar="S", dest="RNAplfold_span_size", type=int, default=150, help=" ")
parser.add_argument("--RNAplfold_unpaired_window", metavar="S", dest="RNAplfold_unpaired_window", type=int, default=8, help=" ")
args = parser.parse_args()

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
config["RNAplfold"]["unpaired_window"]=args.RNAplfold_single_strandedness_cutoff

config["TTS_bed"]["ucsc_dark_gray"]=args.TTS_bed_ucsc_dark_gray
config["TTS_bed"]["dark_gray_stability"]=args.TTS_bed_ucsc_dark_gray_stability

config_yaml = yaml.dump(config)



tmpdir = tempfile.mkdtemp(prefix=args.out_dir + "/" + "3plex_tmp_directory_", dir=".")

with open(tmpdir+"/config.yaml","w") as config_file:
    config_file.write(config_yaml)

bashCommand = "source /etc/profile; conda activate 3plex_v0.1; snakemake --snakefile /opt/3plex/Snakefile -j {jobs} -n {ssRNA}_ssmasked-{dsDNA}.tpx.summary.gz {ssRNA}_ssmasked-{dsDNA}.tpx.stability.gz".format(jobs=args.jobs, ssRNA=args.ssRNA, dsDNA=args.dsDNA)
execute(bashCommand,tmpdir)
