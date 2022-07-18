# 3plex

## Quick start with Docker

Clone the repository to have some testing sequences.
```
mkdir 3plex
wget -O v0.1.2-beta.zip https://github.com/molinerisLab/3plex/archive/refs/tags/v0.1.2-beta.zip
unzip 0.1.2-beta.zip
cd v0.1.2-beta
```

Then run the test using docker pulled from docker hub.
```
docker run -u `id -u` -it --rm -v $PWD:$PWD imolineris/3plex:v0.0.2 $PWD/test/ssRNA.fa $PWD/test/dsDNA.fa $PWD/test_out/
```

Check the output files.
```
ls test_out/*/
```

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

## Running 3plex without docker

The logic of __3plex__ is described in the `docker_context/Snakefile` ad is wrapped by the `3plex.py` script.

## Required input 

- a FASTA file reporting a single sequence (the single stranded RNA sequence). The FASTA header must be the same name of the FASTA file.
- a multi-FASTA file containing one or multiple sequences (the double stranded DNA sequences putatively bound by the RNA via triplex).

In an environment with all the dependencies and scripts available, launch __3plex__ with

```
3plex.py --no_env --snakefile /path/to/3plex/docker_context/Snakefile /path/to/RNA.fa /path/to/DNA.fa /path/to/out_directory
```
