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

—

**Abbreviations:**

`ssRNA`: single-stranded RNA

`dsDNA`: double-stranded DNA

`TFO`: Triplex-Forming Oligo (the ssRNA region that binds the DNA)

`TTS`: Triplex Target Site (the dsDNA region bound by the ssRNA)

`DBD`: DNA Binding Domain (the ssRNA region identified from the overlap of multiple TFOs)


# Web interface

An open beta web interface is available at https://3plex.unito.it/. Please be patient about bugs and report them on [issues](https://github.com/molinerisLab/3plex/issues).


# Command line usage 

3plex can be used to produce raw triplex predictions by downloading and running the 3plex Docker image. For advanced analysis workflows, one can clone the repository and then follow the proposed pipelines listed below. Singularity image is available.

## 1. Run 3plex with Docker

Pull the latest image from the Docker hub:
```
docker pull imolineris/3plex:v1.0
```
How to run:
```
docker run -u `id -u`:`id -g` -it --rm -v $PWD:$PWD imolineris/3plex:v1.0 $PWD/ssRNA.fa $PWD/dsDNA.fa $PWD/results_3plex
```

See the example of [ssRNA.fa](https://github.com/molinerisLab/3plex/blob/main/test/SRA1.fa) and [dsDNA.fa](https://github.com/molinerisLab/3plex/blob/main/test/MANE.GRCh38.v1.1.refseq_genomic.Symbol.tss.1500_500.subset.fa).

###  Required arguments
```
ssRNA.fa          	The RNA sequence in FASTA format. The file must contain only one sequence.
dsDNA.fa          	The DNA sequences in multi-FASTA format.
outdir             	Absolute path to output directory.
```

### Options
 
```
  -j CPUS, --jobs CPUS		Number of parallel threads.
  -l N, --min_length N		Minimum triplex length required. Allowed values: N>=5. 
[ Default: 10 ]
  -e E, --error_rate E  		Maximum percentage of error allowed in a triplex. [ Default: 20 ]
  -s S, --single_strandedness_cutoff S
                    			Percentage of masked RNA nucleotides based on RNAplfold base pairing
probabilities. [ Default: 0 ]
  -c C, --consecutive_errors C
                    			Maximum number of consecutive errors allowed in a triplex. [ Default: 1 ]
  -g G, --guanine_rate G
                    			Minimum percentage of guanines required in a TTS. [ Default: 40 ]
  -r, --filter_repeat   		If enabled, exclude repeat and low complexity regions. [ Default: FALSE ]
  -L M, --max_length M 		Maximum triplex length permitted, M=-1 imply no limits. [ Default: -1 ]
  -t T, --pato_other_parameters T
                    			Additional pato parameters passed as a string (e.g. -t '-mamg 90 -E 4').
3plex output format will not change.
  --pato_simultaneous_sequences SIM_SEQ
                    			Maximum number of sequences that may be processed simultaneously
(less simultaneous sequences equals less memory usage).
[ Default: 256 ]
  --RNAplfold_window_size W_SIZE
                    			RNAplfold: average pair probabilities over windows of specified size.
[ Default: 200 ]
  --RNAplfold_span_size W_SPAN
                    			RNAplfold: maximum separation of a base pair permitted. [ Default: 150 ]
  --RNAplfold_unpaired_window W_LEN
                    			RNAplfold: mean probability that regions of specified W_LEN are unpaired.
[ Default: 8 ]
  --RNA2D_out PATH  		Output the RNAplfold modified z-score in this PATH.
```

> :warning: If you get a *MissingInputException* error make sure your input files are accessible by the user or the group specified in Docker. Pay attention to symlink outside to the mounted volume.

> :warning: If you get no output and docker exit with the 139 status the host linux-kernel version is >= 4.8 and you need to enable vsyscall at startup, see https://helpcenter.onlyoffice.com/installation/mail-enabling-vsyscall.aspx.

### Output

3plex returns one tab-delimited file listing all the TFO:TTS matches: *_tpx.stability.gz_ (see the example of [tpx.stability.gz](https://github.com/molinerisLab/3plex/blob/main/test/results_3plex/SRA1_ssmasked-MANE.GRCh38.v1.1.refseq_genomic.Symbol.tss.1500_500.subset.tpx.stability.gz)).

These are the columns of the file:

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
and a second tab-delimited file reporting a summary triplex score for each dsDNA sequence: *_tpx.summary.add_zeros.gz_ (see the example of [tpx.summary](https://github.com/molinerisLab/3plex/blob/main/test/results_3plex/SRA1_ssmasked-MANE.GRCh38.v1.1.refseq_genomic.Symbol.tss.1500_500.subset.tpx.summary.add_zeros.gz)): 

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
`Stability_best`: the stability score of the most stable TFO:TTS couple among the predicted ones.
`Stability_tot`: the sum of the stability scores resulting from the overlap of the TTSs (highest value).
`Stability_norm`: the dsDNA length normalised Stability_tot score.
`Score_best`: PATO best score (sum of the matches).



## 2. Run 3plex Snakemake workflows

### Dependencies
pato
direnv
viennarna=2.4.7
snakemake=7.8.5
bedtools=2.29.0
gawk

To run 3plex workflows, clone the repository, move inside `3plex`directory and allow direnv:
```
git clone git@github.com:molinerisLab/3plex.git
cd 3plex
direv allow
```

Alternatively, if direnv is not installed, one can manually set the environment variables:
```
export PRJ_ROOT={3plex root directory}
export PATH=$PATH:$PRJ_ROOT/local/bin
```

Create and activate conda environment:
```
conda env create --name 3plex --file=local/envs/3plex.yaml
conda activate 3plex
```

### Raw triplex prediction

This workflow can be used to produce the raw  _tpx.stability.gz_  and _tpx.summary.add_zeros.gz_  files without using Docker.

Move to `dataset/ref_from_sequences` and modify the `from_sequences` section of the `config.yaml` according to your needs following the comments. Then run:

```
snakemake -j N_CORES run_from_sequences
```

### Promoter TPX stability test

This workflow allows the integration of gene expression data to characterise the triplex-forming potential of the investigated ssRNA.

Starting from a list of "universe" genes (e.g., all the expressed genes in the system) and a list of genes of interest (e.g., differentially expressed genes identified upon a lncRNA KD):
retrieve the promoters associated with the genes as annotated in [MANE](http://www.ensembl.org/info/genome/genebuild/mane.html)  
run 3plex with a given ssRNA and the retrieved promoters
compare the stability of the putative triplexes formed with promoters of genes of interest and all the remaining genes with a Mann-Whitney test
perform a [gene set enrichment analysis](https://www.gsea-msigdb.org/gsea/index.jsp) ranking the universe of genes according to their triplex stability score thus computing the significance of the enrichment in promoters with a high or a low TPX stability score and the [leading edge](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading_Edge%20Analysis) table providing a selection of candidate target genes

Move to `dataset/ref_promoter_tpx_stability_test` and modify the `ref_promoter_tpx_stability_test` section in the `config.yaml` according to your needs following the comments. Then run:

```
snakemake -j N_CORES run_promoter_tpx_stability_test
```

Find the results in the `results` directory. To produce the Snakemake report, add the `--report report.html` option when running snakemake.

### Random regions test

Available soon.

## 3. Run 3plex with Singularity

Find it at [3plex/singularity_images/](https://github.com/molinerisLab/3plex/singularity_images/)

—
# Poster

<img src="https://github.com/molinerisLab/3plex/blob/main/.images/poster_cicconetti.png" alt="poster" width="800"/>


