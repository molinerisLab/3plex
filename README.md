<img src="https://github.com/molinerisLab/3plex/blob/main/.images/3plex_logo.png" alt="logo" width="600"/>

# Introduction

3plex is a framework for predicting ssRNA:dsDNA interaction through triple helix (triplex) formation.

3plex integrates the state-of-the-art algorithm for triplex identification with relevant biophysical knowledge on triplex structures: thermal stability derived from triplex denaturation experiments and RNA secondary structure prediction.

The `PATO` algorithm (Amatria-Barral et al., 2023) scans a couple of nucleotide sequences and returns all the triplexes that satisfy a set of user-defined constraints.

3plex enables the exclusion of ssRNA nucleotides from the search based on the secondary structure prediction performed with `RNAplfold` from the ViennaRNA package (Lorentz et al., 2011). Theoretically, a triplex interaction requires a stretch of nucleotides on the ssRNA transcript not to be involved in any further hydrogen bond.

The identified putative triplexes are evaluated according to their thermal stability derived from the `LongTarget` collection of triplex denaturation experiments (He et al., 2015).

Extensive description of the tool can be found in our paper:


>__3plex enables deep computational investigation of triplex forming lncRNAs.__<br>
> Cicconetti C, Lauria A, Proserpio V, Masera M, Tamburrini A, Maldotti M, Oliviero S, Molineris I.<br>
> Comput Struct Biotechnol J. 2023 May 17;21:3091-3102. doi: [10.1016/j.csbj.2023.05.016](https://www.sciencedirect.com/science/article/pii/S2001037023001988). PMID: [37273849](https://pubmed.ncbi.nlm.nih.gov/37273849/); PMCID: PMC10236371.


**Abbreviations**

`ssRNA`: single-stranded RNA

`dsDNA`: double-stranded DNA

`TFO`: Triplex-Forming Oligo (the ssRNA region that binds the DNA)

`TTS`: Triplex Target Site (the dsDNA region bound by the ssRNA)

`DBD`: DNA Binding Domain (the ssRNA region identified from the overlap of multiple TFOs)

---

# Web interface

An open beta web interface is available at https://3plex.unito.it/. Please be patient about bugs and report them on issues.

---

# Run 3plex with Docker

Pull the latest image from the Docker hub:
```
docker pull imolineris/3plex:v1.0
```
How to run:
```
docker run -u `id -u`:`id -g` -it --rm -v $PWD:$PWD imolineris/3plex:v1.0 $PWD/ssRNA.fa $PWD/dsDNA.fa $PWD/results_3plex
```

###  Required arguments
```
ssRNA.fa          	The RNA sequence in FASTA format. The file must contain only one sequence.
dsDNA.fa          	The DNA sequences in multi-FASTA format.
outdir             	Absolute path to output directory.
```

See the example of [ssRNA.fa](https://github.com/molinerisLab/3plex/blob/main/test/SRA1.fa) and [dsDNA.fa](https://github.com/molinerisLab/3plex/blob/main/test/MANE.GRCh38.v1.1.refseq_genomic.Symbol.tss.1500_500.subset.fa).

### Options
 
```
-j CPUS, --jobs CPUS    Number of parallel threads.
-l N, --min_length N    Minimum triplex length required. Allowed values: N>=5. [ Default: 10 ]
-e E, --error_rate E    Maximum percentage of error allowed in a triplex. [ Default: 20 ]
-s S, --single_strandedness_cutoff S
                        Percentage of masked RNA nucleotides based on RNAplfold base pairing probabilities. [ Default: 0 ]
-c C, --consecutive_errors C
                        Maximum number of consecutive errors allowed in a triplex. [ Default: 1 ]
-g G, --guanine_rate G
                        Minimum percentage of guanines required in a TTS. [ Default: 40 ]
-r, --filter_repeat     If enabled, exclude repeat and low complexity regions. [ Default: FALSE ]
-L M, --max_length M    Maximum triplex length permitted, M=-1 imply no limits. [ Default: -1 ]
-t T, --pato_other_parameters T
                        Additional pato parameters passed as a string (e.g. -t '-mamg 90 -E 4'). 3plex output format will not change.
--pato_simultaneous_sequences SIM_SEQ
                        Maximum number of sequences that may be processed simultaneously (less simultaneous sequences equals less memory usage). [ Default: 256 ]
--RNAplfold_window_size W_SIZE
                        RNAplfold: average pair probabilities over windows of specified size. [ Default: 200 ]
--RNAplfold_span_size W_SPAN
                        RNAplfold: maximum separation of a base pair permitted. [ Default: 150 ]
--RNAplfold_unpaired_window W_LEN
                        RNAplfold: mean probability that regions of specified W_LEN are unpaired. [ Default: 8 ]
--RNA2D_out PATH        Output the RNAplfold modified z-score in this PATH.
```

> :warning: If you get a *MissingInputException* error make sure your input files are accessible by the user or the group specified in Docker. Pay attention to symlink outside to the mounted volume.

> :warning: If you get no output and docker exit with the 139 status the host linux-kernel version is >= 4.8 and you need to enable vsyscall at startup, see https://helpcenter.onlyoffice.com/installation/mail-enabling-vsyscall.aspx.

### Output

3plex returns one tab-delimited file listing all the TFO:TTS matches: *_tpx.stability.gz_ (see the example of [tpx.stability](https://github.com/molinerisLab/3plex/blob/main/test/results_3plex/SRA1_ssmasked-MANE.GRCh38.v1.1.refseq_genomic.Symbol.tss.1500_500.subset.tpx.stability.txt)).

```
1   	Sequence_ID
2   	TFO_start
3   	TFO_end
4   	Duplex_ID
5   	TTS_start
6   	TTS_end
7   	Score
8   	Error_rate
9   	Errors
10  	Motif
11  	Strand
12  	Orientation
13  	Guanine_rate
14  	Stability
15  	aln1
16  	aln2
17  	aln3
18  	aln4
```

and a second tab-delimited file reporting a summary triplex score for each dsDNA sequence: *_tpx.summary.add_zeros.gz_ (see the example of [tpx.summary](https://github.com/molinerisLab/3plex/blob/main/test/results_3plex/SRA1_ssmasked-MANE.GRCh38.v1.1.refseq_genomic.Symbol.tss.1500_500.subset.tpx.summary.add_zeros.txt)): 

```
1   	Duplex_ID
2   	Sequence_ID
3   	Total(abs)
4   	Total(rel)
5   	GA(abs)
6   	GA(rel)
7   	TC(abs)
8   	TC(rel)
9   	GT(abs)
10  	GT(rel)
11  	Duplex_length
12  	Stability_best
13  	Stability_tot
14  	Score_best
15  	Stability_norm
```

**Available triplex scores**

`Stability_best`: the stability score of the most stable TFO:TTS couple among the predicted ones.

`Stability_tot`: the sum of the stability scores of overlapping TTSs (highest value).

`Stability_norm`: the dsDNA length normalised Stability_tot score.

`Score_best`: PATO best score (sum of the matches).

---

# Run 3plex with Singularity

Find the Singularity image at [3plex/singularity_images/](https://github.com/molinerisLab/3plex/singularity_images/)

---

# 3plex Snakemake workflows

The following section illustrates how to perform 3plex downstream analyses to additionally evaluate the significance of the triplex-forming capability of an ssRNA.

The workflows are handled with [Snakemake](https://snakemake.github.io/) and structured with [dap](https://github.com/molinerisLab/dap).

## Dependencies

```
pato
dap
gawk
```
Find here the [installation guide for PATO](https://github.com/UDC-GAC/pato).

## Setup

We suggest to install [direnv](https://direnv.net/). Alternatively, one can manually set the environment variables:
```
export PRJ_ROOT={3plex_root_directory}
export PATH=$PATH:$PRJ_ROOT/local/bin
```

To run the workflows clone the repository, move inside `3plex` directory and, if installed, allow direnv:
```
git clone git@github.com:molinerisLab/3plex.git
cd 3plex
direnv allow
```

Create and activate the conda environment:
```
conda env create --file=local/envs/3plex.yaml
conda activate 3plex_Env
```
Then move to the directory that contains the first version of your analysis:
```
cd v1
```
To run a different version of your analysis you can clone v1 into v2 by running `dap clone v1 v2`.

## Raw triplex prediction

This workflow produces the raw  _tpx.stability.gz_ and _tpx.summary.add_zeros.gz_ files previously described. An ssRNA FASTA file and a dsDNA multi-FASTA or BED file are required.

Specify the `ssRNA` and `dsDNA` paths in the `config.yaml`. Then run:

```
snakemake -j CORES run_raw_tpx_prediction
```

## Promoter TPX stability test

This workflow allows the integration of gene expression data to characterize the triplex-forming potential of the investigated ssRNA.

Starting from a list of the "universe of genes" (e.g., all the expressed genes in the system) and a list of genes of interest (e.g., differentially expressed genes identified upon an ssRNA KD):

1. retrieve the promoters associated with the universe of genes (we suggest to refer to [MANE](http://www.ensembl.org/info/genome/genebuild/mane.html), but a custom selection of promoters is allowed);
2. run 3plex to find the putative triplexes formed by the ssRNA and the set of promoters;
3. compare the stability of the putative triplexes formed with promoters of genes of interest with all the remaining genes (Wilcoxon's test);
4. perform a [gene set enrichment analysis](https://www.gsea-msigdb.org/gsea/index.jsp) ranking the universe of genes according to their triplex stability score thus computing the significance of the enrichment in promoters with a high or a low stability score. The [leading edge](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading_Edge%20Analysis) table provides a selection of candidate target genes.

Modify the `Promoter stability test` section in the `config.yaml` according to your needs following the comments. Then run:

```
snakemake -j CORES run_promoter_tpx_stability_test
```

Find the results in the `results` directory. 

To produce an HTML report, run the following command after the workflow is completed:

```
snakemake -j CORES run_promoter_tpx_stability_test --report report.html
```

## Random regions test

This workflow tests the triplex-forming capability of each ssRNA portion, namely the DBDs. 

This result is achieved by comparing the stability of the DBDs' putative triplexes identified considering a set of target regions (e.g., ChIRP-seq data) with a null distribution built on the stability scores obtained testing randomized genomic regions N times.

Operatively:

1. run 3plex on the set of target regions and define the DBDs as the union of the overlapping TFOs;
2. generate N times the same number of regions with bedtools shuffle, excluding blacklist regions. 3plex is run each time to retrieve the TTS count and stability of the triplexes;
3. Compute the empirical p-value for each DBD by counting the number of times the upper quartile of the DBD's putative triplexes' score is lower than the upper quartile obtained with random genomic regions.

Modify the `Random region test` section in the `config.yaml` according to your needs following the comments. Then run:

```
snakemake -j CORES run_random_region_test
```

Find the results in the `results` directory. 

---

# Poster

<img src="https://github.com/molinerisLab/3plex/blob/main/.images/poster_cicconetti.png" alt="poster" width="800"/>
