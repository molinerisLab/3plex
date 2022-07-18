# 3plex

## Introduction

3plex is a software that predict the interaction between a single strand RNA (ssRNA) with a double strand DNA region (dsDNA) through triple helix formation (ssRNA:dsDNA TPX). 

`Triplexator` algorithm (Buske et al., 2011)  is used to scan a couple of input nucleotide sequences and to return all the TPX that satisfy a set of user-defined constraints. We lowered the minimum length required for the definition of a TPX from 10 to 5 ([see our paper results](https://doi.org/10.1101/2022.07.06.496678)).

The identified putative TPX are scored according to their thermal stability derived from `LongTarget` collection of thermal denaturation experiments (He et al., 2015).

3plex integrates `RNAplfold` from the ViennaRNA package (Lorentz et al., 2011) to consider the RNA secondary structure information in the definition of a TPX.

---

## Quick start with Docker

By downloading the software release you save the testing sequences.
```
mkdir 3plex
wget -O v0.1.2-beta.zip https://github.com/molinerisLab/3plex/archive/refs/tags/v0.1.2-beta.zip
unzip v0.1.2-beta.zip
cd v0.1.2-beta
```

Then run the test using docker pulled from docker hub.
```
docker run -u `id -u` -it --rm -v $PWD:$PWD imolineris/3plex:v0.1.2-beta $PWD/test/ssRNA.fa $PWD/test/dsDNA.fa $PWD/test_out/
```

Check the output files.
```
ls test_out/*/
```

---

## How to install

To simplify the installation, we advise you to use the docker image. See the dockerfile if you need to modify the software.

### Dependencies

- triplexator ()
- viennarna=2.4.7
- snakemake=7.8.5
- bedtools=2.29.0
- gawk

### Executable scripts

Clone the repository
```git clone```

then add the scripts contained in the `docker_context` directory to PATH variable
```
cd 3plex
export PATH:$PATH:$PWD/docker_context
```

---

## Running 3plex without docker

The logic of __3plex__ is described in the `docker_context/Snakefile` ad is wrapped by the `3plex.py` script.

In an environment with all the dependencies and scripts available, launch __3plex__ with

```
3plex.py --no_env --snakefile /path/to/3plex/docker_context/Snakefile /path/to/RNA.fa /path/to/DNA.fa /path/to/out_directory
```

---

## 3plex usage

### Required input 

- a FASTA file reporting the single ssRNA sequence ([Example of ssRNA.fa](https://github.com/molinerisLab/3plex/blob/main/test/ssRNA.fa))
- a multi-FASTA file containing one or multiple dsDNA sequences ([Example of dsDNA.fa](https://github.com/molinerisLab/3plex/blob/main/test/dsDNA.fa)).
- a defined output directory


### Output format

```
```


### Options
 
```
  -h, --help            show this help message and exit
  --bed dsDNA.bed       Genomic coordiantes of the DNA sequences in bed
                        format, the 4th column must contain the same
                        identifiers of sequences in dsDNA.fa (default: None)
  -j CPUS, --jobs CPUS  Number of parallel threads. (default: 1)
  -l N, --min_length N  Minimum triplex length required. Allowed values: N>5.
                        (default: 10)
  -e E, --error_rate E  Maximal percentage of error allowed in a triplex.
                        (default: 20)
  -s S, --single_strandedness_cutoff S
                        Percentage of masked RNA nucleotides based on
                        RNAplfold base pairing probabilities. (default: 0)
  -c C, --consecutive_errors C
                        Maximum number of consecutive errors allowed in a
                        triplex. (default: 3)
  -g G, --guanine_rate G
                        Minimal guanine percentage required in any TTS.
                        (default: 10)
  -r R, --filter_repeat R
                        If enabled, exclude repeat and low complexity regions.
                        (default: False)
  -L M, --max_length M  Maximum triplex length permitted, M=-1 imply no
                        limits. (default: -1)
  -t T, --triplexator_other_parameters T
                        Additional triplexator parameters passed as a sting
                        (e.g. -t '-mamg 90 -E 4'). Triplexator output format
                        will not change. (default: )
  --ucsc_dark_gray G    TTS bed UCSC dark gray (default: 843)
  --dark_gray_stability G
                        10% of TTS in paper. (default: 43)
  --RNAplfold_window_size S
                        RNAplfold: average pair probabilities over windows of
                        specified size. (default: 200)
  --RNAplfold_span_size S
                        RNAplfold: maximum separation of a base pair
                        permitted. (default: 150)
  --RNAplfold_unpaired_window S
                        RNAplfold: mean probability that regions of specified
                        length are unpaired. (default: 8)
  --snakefile file      (default: /opt/3plex/Snakefile)
  --no_env NO_ENV       Do not load the conda environment, useful when running
                        the script outside of the docker image. (default:
                        False)
```
