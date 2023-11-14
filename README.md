<img src="https://github.com/molinerisLab/3plex/blob/main/www/3plex_logo.png" alt="logo" width="600"/>

# Introduction

3plex is a framework for predicting ssRNA:dsDNA interaction through triple helix (triplex) formation.

3plex integrates the state-of-the-art algorithm for triplex identification with relevant biophysical knowledge on triplex structures: thermal stability derived from triplex denaturation experiments and RNA secondary structure prediction.

`PATO` algorithm (Amatria-Barral et al., 2023) scans a couple of nucleotide sequences and returns all the triplex that satisfy a set of user-defined constraints.

The identified putative triplexes are scored according to their thermal stability derived from the `LongTarget` collection of triplex denaturation experiments (He et al., 2015).

3plex integrates `RNAplfold` from the ViennaRNA package (Lorentz et al., 2011) to enable the exclusion of the ssRNA regions not available to the binding from the prediction.

Extensive description of the tool can be found in our paper:


>__3plex enables deep computational investigation of triplex forming lncRNAs.__<br>
> Cicconetti C, Lauria A, Proserpio V, Masera M, Tamburrini A, Maldotti M, Oliviero S, Molineris I.<br>
> Comput Struct Biotechnol J. 2023 May 17;21:3091-3102. doi: [10.1016/j.csbj.2023.05.016](https://www.sciencedirect.com/science/article/pii/S2001037023001988). PMID: [37273849](https://pubmed.ncbi.nlm.nih.gov/37273849/); PMCID: PMC10236371.

---

# Web interface

A web interface to easily work with 3plex is coming soon!

---

# Quick start with Docker

Just to have some testing sequences, clone the repository:
```
mkdir 3plex
wget -O v0.1.2-beta.zip https://github.com/molinerisLab/3plex/archive/refs/tags/v0.1.2-beta.zip
unzip v0.1.2-beta.zip
cd v0.1.2-beta
```

You can test the tool using:

 * `test/ssRNA.fa`
 * `test/dsDNA.fa`

You do not need to use directly the code in the repository or install dependencies. You can run a test using the image pulled from docker hub.
```
cd 3plex;
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

> :warning: If you get no output and docker exit with the 139 status the host linux-kernel version is >= 4.8 and you need to enable vsyscall at startup, see https://helpcenter.onlyoffice.com/installation/mail-enabling-vsyscall.aspx.
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
  --bed dsDNA.bed       Genomic coordinates of the DNA sequences in bed
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
                        If enabled, exclude repeat and low-complexity regions.
                        (default: False)
  -L M, --max_length M  Maximum triplex length permitted, M=-1 imply no
                        limits. (default: -1)
  -t T, --triplexator_other_parameters T
                        Additional triplexator parameters passed as a string
                        (e.g. -t '-mamg 90 -E 4'). Triplexator output format
                        will not change. (default: )
  --ucsc_dark_gray G    TTS bed UCSC dark grey (default: 843)
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
# Singularity images

Find it at [3plex/singularity_images/](https://github.com/molinerisLab/3plex/singularity_images/)

---
## Poster

<img src="https://github.com/molinerisLab/3plex/blob/main/www/poster_cicconetti.png" alt="poster" width="1000"/>
