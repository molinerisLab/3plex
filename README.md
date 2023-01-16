<img src="https://github.com/molinerisLab/3plex/blob/main/www/3plex_logo.png" alt="logo" width="600"/>

# Introduction

3plex is a software that predict the interaction between a single strand RNA (ssRNA) with a double strand DNA region (dsDNA) through triple helix formation (ssRNA:dsDNA TPX). 

`Triplexator` algorithm (Buske et al., 2011)  is used to scan a couple of input nucleotide sequences and to return all the TPX that satisfy a set of user-defined constraints. We lowered the minimum length required for the definition of a TPX from 10 to 5 ([see our paper results](https://doi.org/10.1101/2022.07.06.496678)).

The identified putative TPX are scored according to their thermal stability derived from `LongTarget` collection of thermal denaturation experiments (He et al., 2015).

3plex integrates `RNAplfold` from the ViennaRNA package (Lorentz et al., 2011) to consider the RNA secondary structure information in the definition of a TPX.

Extensive description of the tool can be foud in out paper:


> __RNABSdb and 3plex enable deep computational investigation of triplex forming lncRNAs__</br>
> Chiara Cicconetti, Andrea Lauria, Valentina Proserpio, Annalaura Tamburrini, Mara Maldotti, Salvatore Oliviero, Ivan Molineris</br>
> bioRxiv 2022.07.06.496678; doi: https://doi.org/10.1101/2022.07.06.496678


---

# Quick start with Docker

Clone the repository.
```
mkdir 3plex
wget -O v0.1.2-beta.zip https://github.com/molinerisLab/3plex/archive/refs/tags/v0.1.2-beta.zip
unzip v0.1.2-beta.zip
cd v0.1.2-beta
```
This is just to have some testing sequences:

 * `test/ssRNA.fa`
 * `test/dsDNA.fa`

You do not need to use directly the code in the repository or install dependencies. You can run a test using the image pulled from docker hub.
```
docker run -u `id -u`:`id -g` -it --rm -v $PWD:$PWD imolineris/3plex:v0.1.2-beta $PWD/test/ssRNA.fa $PWD/test/dsDNA.fa $PWD/test_out/
```

Check the output files.
```
ls test_out/*/
```

To see the option list:
```
docker run -u `id -u`:`id -g` -it --rm -v $PWD:$PWD imolineris/3plex:v0.1.2-beta -h
```

To run 3plex on your data just change the test fasta files with the ones you are interested in.

> :warning: If you get a *MissingInputException* error make sure your input files are accessible by the user or the group specified in Docker.

---

# 3plex usage

## Required inputs 

- a FASTA file reporting the single ssRNA sequence ([Example of ssRNA.fa](https://github.com/molinerisLab/3plex/blob/main/test/ssRNA.fa))
- a multi-FASTA file containing one or multiple dsDNA sequences ([Example of dsDNA.fa](https://github.com/molinerisLab/3plex/blob/main/test/dsDNA.fa)).
- a defined output directory


## Outputs format

- a tab delimited file, named _tpx.stability_ ([Example of tpx.stability](https://github.com/molinerisLab/3plex/blob/main/test/test_out/Z85996.1_ssmasked-dsDNA.tpx.stability)) reporting:

```
1       Sequence_ID
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
15      aln1
16      aln2
17      aln3
18      aln4
```

- a tab delimited file, named _tpx.summary_ ([Example of tpx.summary](https://github.com/molinerisLab/3plex/blob/main/test/test_out/Z85996.1_ssmasked-dsDNA.tpx.summary)) reporting:

```
1       Duplex_ID
2       Sequence_ID
3       Total(abs)
4       Total(rel)
5       GA(abs)
6       GA(rel)
7       TC(abs)
8       TC(rel)
9       GT(abs)
10      GT(rel)
11      Duplex_length
12      Stability_best
13      Stability_tot
14      Score_best
15      Stability_norm
```


## Options
 
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

---


# How to install

To simplify the installation, we advise you to use the docker image. See the dockerfile if you need to modify the software.

## Dependencies

- triplexator
- viennarna=2.4.7
- snakemake=7.8.5
- bedtools=2.29.0
- gawk

## Executable scripts

Clone the repository
```git clone```

then add the scripts contained in the `docker_context` directory to PATH variable
```
cd 3plex
export PATH:$PATH:$PWD/docker_context
```

---

# Running 3plex without docker

The logic of __3plex__ is described in the `docker_context/Snakefile` ad is wrapped by the `3plex.py` script.

In an environment with all the dependencies and scripts available, launch __3plex__ with

```
3plex.py --no_env --snakefile /path/to/3plex/docker_context/Snakefile /path/to/RNA.fa /path/to/DNA.fa /path/to/out_directory
```

---
# Singularity images

Find it at [3plex/singularity_images/](https://github.com/molinerisLab/3plex/singularity_images/)
