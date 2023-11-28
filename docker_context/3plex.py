#!/usr/bin/env python3
import os
import re
import argparse
import yaml
import tempfile
import subprocess

REMOVE_TMP_DIR = True
PRINT_ACTIONS = True

def print_log(text):
    if (PRINT_ACTIONS):
        print(text)

def execute(cmd, cwd):
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
    # parser
    parser = argparse.ArgumentParser(
    description='''
Given a ssRNA sequence and a list of dsDNA sequences, 3plex finds all the triplexes that satisfy the constraints and computes a thermal stability score for the interaction.
The RNA secondary structure prediction can be used to exclude RNA nucleotides from the search.
    ''',
    formatter_class= argparse.RawTextHelpFormatter,
    epilog='''
Please cite:
    Cicconetti C, Lauria A, Proserpio V, et al. 3plex enables deep computational investigation of triplex forming lncRNAs. Comput Struct Biotechnol J. 2023;21:3091-3102.
    
Have a good 3plex journey!

    '''
    )
    # positional arguments
    parser.add_argument("ssRNA", metavar="ssRNA.fa", 
                        help="The RNA sequence in FASTA format. The file must contain only one sequence.")
    parser.add_argument("dsDNA", metavar="dsDNA.fa", 
                        help="The DNA sequences in multi-FASTA format.")
    parser.add_argument("out_dir", metavar="PATH",
                        help="Absolute path to output directory.")

    # optional arguments
    parser.add_argument("-j", "--jobs", metavar="CPUS", dest="jobs", type=int, default=1, 
                        help="Number of parallel threads.")
    parser.add_argument("-l", "--min_length", metavar="N", dest="pato_min_length", type=minlen_type, default=10,          
                        help="Minimum triplex length required. Allowed values: N>=5. [ Default: 10 ]")
    parser.add_argument("-e", "--error_rate", metavar="E", dest="pato_error_rate", type=int, default=20, choices=range(21), 
                        help="Maximum percentage of error allowed in a triplex. [ Default: 20 ]")
    parser.add_argument("-s", "--single_strandedness_cutoff", metavar="S", dest="RNAplfold_single_strandedness_cutoff", type=int, default=0, choices=range(100), 
                        help="Percentage of masked RNA nucleotides based on RNAplfold base pairing probabilities. [ Default: 0 ]")
    parser.add_argument("-c", "--consecutive_errors", metavar="C", dest="pato_consecutive_errors", type=int, default=1, 
                        help="Maximum number of consecutive errors allowed in a triplex. [ Default: 1 ]")
    parser.add_argument("-g", "--guanine_rate", metavar="G", dest="pato_guanine_rate", type=int, default=40, choices=range(100), 
                        help="Minimum percentage of guanines required in a TTS. [ Default: 40 ]")
    parser.add_argument("-r", "--filter_repeat", dest="pato_filter_repeat", action="store_true",
                        help="If enabled, exclude repeat and low complexity regions. [ Default: FALSE ]")
    parser.add_argument("-L", "--max_length", metavar="M", dest="pato_max_length", type=int, default=-1, 
                        help="Maximum triplex length permitted, M=-1 imply no limits. [ Default: -1 ]")
    parser.add_argument("-t", "--pato_other_parameters", metavar="T", dest="pato_other_parameters", type=str, default="", 
                        help="Additional pato parameters passed as a string (e.g. -t '-mamg 90 -E 4'). 3plex output format will not change.")
    parser.add_argument("--pato_simultaneous_sequences", metavar="SIM_SEQ", type=int, default=256,
                        help="Maximum number of sequences that may be processed simultaneously (less simultaneous sequences equals less memory usage). [ Default: 256 ]"),
    parser.add_argument("--RNAplfold_window_size", metavar="W_SIZE", dest="RNAplfold_window_size", type=int, default=200, 
                        help="RNAplfold: average pair probabilities over windows of specified size. [ Default: 200 ]")
    parser.add_argument("--RNAplfold_span_size", metavar="W_SPAN", dest="RNAplfold_span_size", type=int, default=150, 
                        help="RNAplfold: maximum separation of a base pair permitted. [ Default: 150 ]")
    parser.add_argument("--RNAplfold_unpaired_window", metavar="W_LEN", dest="RNAplfold_unpaired_window", type=int, default=8, 
                        help="RNAplfold: mean probability that regions of specified W_LEN are unpaired. [ Default: 8 ]")
    parser.add_argument("--RNA2D_out", metavar="PATH", dest="RNA2D_out", type=str, default=None, 
                        help="Output the RNAplfold modified z-score in this PATH.")
    args = parser.parse_args()

    if args.pato_filter_repeat:
        args.pato_filter_repeat="on"
    else:
        args.pato_filter_repeat="off"

    # build config object
    config={}

    config["pato"]={}
    config["RNAplfold"]={}

    config["pato"]["min_length"]=args.pato_min_length
    config["pato"]["max_length"]=args.pato_max_length
    config["pato"]["error_rate"]=args.pato_error_rate
    config["pato"]["guanine_rate"]=args.pato_guanine_rate
    config["pato"]["filter_repeat"]=args.pato_filter_repeat
    config["pato"]["consecutive_errors"]=args.pato_consecutive_errors
    config["pato"]["other_parameters"]=args.pato_other_parameters
    config["pato"]["sim_sequences"]=args.pato_simultaneous_sequences

    config["RNAplfold"]["window_size"]=args.RNAplfold_window_size
    config["RNAplfold"]["span_size"]=args.RNAplfold_span_size
    config["RNAplfold"]["unpaired_window"]=args.RNAplfold_unpaired_window
    config["RNAplfold"]["single_strandedness_cutoff"]=args.RNAplfold_single_strandedness_cutoff

    # create output dir
    print_log("Preparing execution environment...")
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    tmpdir = tempfile.mkdtemp(prefix=args.out_dir + "/" + "3plex_tmp_directory_", dir=".")

    # export config
    config_yaml = yaml.dump(config)
    with open(tmpdir+"/config.yaml","w") as config_file:	
        config_file.write(config_yaml)

    # prepare directory
    # link to ssRNA and sdDNA files
    ssRNA_name=get_ssRNA_seq_name(args.ssRNA)
    dsDNA_name = os.path.splitext(os.path.basename(args.dsDNA))[0]
    os.symlink(os.path.abspath(args.ssRNA), f"{tmpdir}/{ssRNA_name}.fa")
    os.symlink(os.path.abspath(args.dsDNA), f"{tmpdir}/{dsDNA_name}.fa")
    # link to Snakefile
    os.symlink(os.path.join(os.getcwd(), "Snakefile"), os.path.join(tmpdir, "Snakefile"))
    print_log(f"Input files:\n>{args.ssRNA}\n>{args.dsDNA}")

    # prepare bash command
    bashCommand = """
    export MAMBA_ROOT_PREFIX=/3plex/mamba;
    eval "$(micromamba shell hook --shell bash)" && micromamba activate && micromamba activate 3plex;
    export PATH=$PATH:/3plex/bin
    """ 
    bashCommand+=f"""
    cd {tmpdir};
snakemake -p -c{args.jobs} \
    {ssRNA_name}_ssmasked-{dsDNA_name}.tpx.summary.add_zeros.gz \
    {ssRNA_name}_ssmasked-{dsDNA_name}.tpx.stability.gz 2>>{tmpdir}/tpx.log;
    mv {ssRNA_name}_ssmasked-{dsDNA_name}.tpx.summary.add_zeros.gz {ssRNA_name}_ssmasked-{dsDNA_name}.tpx.stability.gz ../;
    mv tpx.log ../;
    mv {ssRNA_name}_ssmasked.fa ../;
    mv config.yaml ../;
"""

    if args.RNA2D_out is not None:
        bashCommand+=f"\n mv RNAplfold/{ssRNA_name}_lunp.modif_zscore ../{args.RNA2D_out};"
    
    # remove tmp dir
    if (REMOVE_TMP_DIR):
        bashCommand+=f"""
        rm -rf {tmpdir}"""
    
    # execute command
    print_log("Running 3plex...")
    execute(bashCommand, tmpdir)
    print_log(f"Done. Output files' path: {args.out_dir}")

if __name__ == "__main__":
    main()