# 3plex

## Qick start with Docker

Download the software release to have the testing sequences
```
mkdir 3plex
wget -O v0.1.1-beta.zip https://github.com/molinerisLab/3plex/archive/refs/tags/v0.1.1-beta.zip
unzip 0.1.1-beta.zip
cd v0.1.1-beta
```

Then run the test using docker pulled from docher hub.
```
docker run -u `id -u` -it --rm -v $PWD:$PWD imolineris/3plex:v0.0.2 $PWD/test/ssRNA.fa $PWD/test/dsDNA.fa $PWD/test_out/
```

Check the outuput files

```
ls test_out/*/
```

## How to install

To simplyfi the installation we advise to use the docker image, see the dockerfile if you need to modify the software.

### Dependencies

- triplexator ()
- viennarna=2.4.7
- snakemake=7.8.5
- bedtools=2.29.0
- gawk

### Execuitable scripts

Clone the repository
```git clone```
and put the scripts in the directory `docker_context` in the PATH
```
cd 3plex
export PATH:$PATH:$PWD/docker_context
```

## Running 3plex without docker

The logic of __3plex__ is described in the `docker_context/Snakefile` ad is wrapped by the `3plex.py` script.

Required input 

- a file containing a single sequence, meant to be a RNA sequence, in the following is the name of this file,
- a file containing one or many sequences, meant to be DNA sequences putatively bound by the RNA via triple helices.

In an environemnt with all the dependencies and scripts available, lounch __3plex__ with

```
3plex.py --no_env --snakefile /path/to/3plex/docker_context/Snakefile /path/to/RNA.fa /path/to/DNA.fa /path/to/out_directory
```
