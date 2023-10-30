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

def minlen_type(x):
    x = int(x)
    if x < 5:
        raise argparse.ArgumentTypeError("Minimum length allowed is 5")
    return x

def get_ssRNA_seq_name(ssRNA):
	with open(ssRNA) as fasta:
		header=fasta.readline()
		m=re.match(r'>([^\s]+)',header)
		return m.group(1)


def main():
    #Parse arguments
    parser = argparse.ArgumentParser(
    description="Given an RNA sequence and a list of DNA sequences, 3plex finds all the triplexes that satisfy the constraints and computes a thermal stability score for the interaction. The RNA secondary structure prediction can be used to exclude RNA nucleotides from the search.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    epilog="Please cite: Cicconetti C, Lauria A, Proserpio V, et al. 3plex enables deep computational investigation of triplex forming lncRNAs. Comput Struct Biotechnol J. 2023;21:3091-3102.\n\nHave a good 3plex journey!")

    # positional
    parser.add_argument("ssRNA", metavar="ssRNA.fa", 
                        help="The RNA sequence in FASTA format. The file must contain only one sequence.")
    parser.add_argument("dsDNA", metavar="dsDNA.fa", 
                        help="The DNA sequences in multi-FASTA format.")
    parser.add_argument("out_dir", metavar="out_dir", 
                        help="Relative path to output directory.")

    # optional
    parser.add_argument("--bed", metavar="dsDNA.bed", dest="dsDNA_bed", type=str, default=None, 
                        help="Genomic coordinates of the DNA sequences in BED format. The 4th column must contain the same identifiers of the sequences in dsDNA.fa")
    parser.add_argument("-j", "--jobs", metavar="CPUS", dest="jobs", type=int, default=1, 
                        help="Number of parallel threads.")
    parser.add_argument("-l", "--min_length", metavar="N", dest="triplexator_min_length", type=minlen_type, default=8,          
                        help="Minimum triplex length required. Allowed values: N>=5.")
    parser.add_argument("-e", "--error_rate", metavar="E", dest="triplexator_error_rate", type=int, default=20, choices=range(21), 
                        help="Maximum percentage of error allowed in a triplex.")
    parser.add_argument("-s", "--single_strandedness_cutoff", metavar="S", dest="RNAplfold_single_strandedness_cutoff", type=int, default=0, choices=range(100), 
                        help="Percentage of masked RNA nucleotides based on RNAplfold base pairing probabilities.")
    parser.add_argument("-c", "--consecutive_errors", metavar="C", dest="triplexator_consecutive_errors", type=int, default=1, 
                        help="Maximum number of consecutive errors allowed in a triplex.")
    parser.add_argument("-g", "--guanine_rate", metavar="G", dest="triplexator_guanine_rate", type=int, default=40, choices=range(100), 
                        help="Minimum percentage of guanines required in a TTS.")
    parser.add_argument("-r", "--filter_repeat", dest="triplexator_filter_repeat", action="store_true",
                        help="If enabled, exclude repeat and low complexity regions.")
    parser.add_argument("-L", "--max_length", metavar="M", dest="triplexator_max_length", type=int, default=-1, 
                        help="Maximum triplex length permitted, M=-1 imply no limits.")
    parser.add_argument("-t", "--triplexator_other_parameters", metavar="T", dest="triplexator_other_parameters", type=str, default="", 
                        help="Additional triplexator parameters passed as a sting (e.g. -t '-mamg 90 -E 4'). 3plex output format will not change.")
    parser.add_argument("--pato_simultaneus_sequences", type=int, default=256,
                        help="Maximum number of sequences that may be processed simultaneously (less simultaneous sequences equals less memory usage)."),
    parser.add_argument("--ucsc_dark_gray", metavar="G", dest="TTS_bed_ucsc_dark_gray", type=int, default=843, choices=range(1000), 
                        help="TTS BED UCSC dark gray.")
    parser.add_argument("--dark_gray_stability", metavar="G", dest="TTS_bed_ucsc_dark_gray_stability", type=int, default=43, 
                        help="10%% of TTS in paper.")
    parser.add_argument("--RNAplfold_window_size", metavar="S", dest="RNAplfold_window_size", type=int, default=200, 
                        help="RNAplfold: average pair probabilities over windows of specified size.")
    parser.add_argument("--RNAplfold_span_size", metavar="S", dest="RNAplfold_span_size", type=int, default=150, 
                        help="RNAplfold: maximum separation of a base pair permitted.")
    parser.add_argument("--RNAplfold_unpaired_window", metavar="S", dest="RNAplfold_unpaired_window", type=int, default=8, 
                        help="RNAplfold: mean probability that regions of specified length are unpaired.")
    parser.add_argument("--RNA2D_out", metavar="PATH", dest="RNA2D_out", type=str, default=None, 
                        help="Output the RNAplfold modified z-score in this PATH.")
    parser.add_argument("--snakefile", metavar="file", dest="snakefile", type=str, default="/opt/3plex/Snakefile", 
                        help="Override the default Snakefile using the one specified.")
    args = parser.parse_args()

    if args.triplexator_filter_repeat:
        args.triplexator_filter_repeat="on"
    else:
        args.triplexator_filter_repeat="off"

    #Build config object
    config={}
    config["triplexator"]={}
    config["RNAplfold"]={}
    config["TTS_bed"]={}
    config["triplexator_other"]={}

    config["triplexator"]["min_length"]=args.triplexator_min_length
    config["triplexator"]["max_length"]=args.triplexator_max_length
    config["triplexator"]["error_rate"]=args.triplexator_error_rate
    config["triplexator"]["guanine_rate"]=args.triplexator_guanine_rate
    config["triplexator"]["filter_repeat"]=args.triplexator_filter_repeat
    config["triplexator"]["consecutive_errors"]=args.triplexator_consecutive_errors

    config["triplexator_other"]["other_parameters"]=args.triplexator_other_parameters

    config["RNAplfold"]["window_size"]=args.RNAplfold_window_size
    config["RNAplfold"]["span_size"]=args.RNAplfold_span_size
    config["RNAplfold"]["unpaired_window"]=args.RNAplfold_unpaired_window
    config["RNAplfold"]["single_strandedness_cutoff"]=args.RNAplfold_single_strandedness_cutoff

    config["TTS_bed"]["ucsc_dark_gray"]=args.TTS_bed_ucsc_dark_gray
    config["TTS_bed"]["dark_gray_stability"]=args.TTS_bed_ucsc_dark_gray_stability

    config["pato_simultaneus_sequences"]= args.pato_simultaneus_sequences

    #Create output dir
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    tmpdir = tempfile.mkdtemp(prefix=args.out_dir + "/" + "3plex_tmp_directory_", dir=".")
    tmpdir = os.path.abspath(tmpdir)

    #Esporta config
    config_yaml = yaml.dump(config)
    with open(tmpdir+"/config.yaml","w") as config_file:	
        config_file.write(config_yaml)

    #Prepare directory:
    #Link ssRNA and sdDNA files
    ssRNA_name=get_ssRNA_seq_name(args.ssRNA)
    os.symlink(os.path.abspath(args.ssRNA), f"{tmpdir}/{ssRNA_name}.fa")
    os.symlink(os.path.abspath(args.dsDNA), os.path.join(tmpdir, os.path.basename(args.dsDNA)))
    dsDNA_name = os.path.splitext(os.path.basename(args.dsDNA))[0]
    #Link snakefile
    os.symlink(os.path.join(os. getcwd(), "Snakefile"), os.path.join(tmpdir, "Snakefile"))
    #TODO test that it works with dsDNA bed + why dsDNA as fasta is needed anyway?
    if args.dsDNA_bed:
        os.symlink(args.dsDNA_bed, os.path.join(tmpdir, os.path.basename(dsDNA_bed)))

    #Prepare bash command
    bashCommand = """
    eval "$(micromamba shell hook --shell bash)" && source ~/.bashrc && micromamba activate;
    micromamba activate 3plex;
    export PATH=$PATH:/3plex/bin
    """ 
    bashCommand+=f"""
    cd {os.path.abspath(tmpdir)};
snakemake -c{args.jobs} \
    {ssRNA_name}_ssmasked-{dsDNA_name}.tpx.summary.add_zeros.gz \
    {ssRNA_name}_ssmasked-{dsDNA_name}.tpx.stability.gz >> {tmpdir}/STDOUT 2>>{tmpdir}/STDERR;
    mv {ssRNA_name}_ssmasked-{dsDNA_name}.tpx.summary.add_zeros.gz {ssRNA_name}_ssmasked-{dsDNA_name}.tpx.stability.gz ../;
    mv STDOUT STDERR ../;
"""

    if args.RNA2D_out is not None:
        bashCommand+=f"\n mv RNAplfold/{ssRNA_name}_lunp.unpairedWindow.modif_zscore ../{args.RNA2D_out};"
    
    bashCommand+=f"""
    rm -rf {tmpdir}"""
    
    #Execute command
    execute(bashCommand, tmpdir)

if __name__ == "__main__":
    main()