# 3plex

## Docker usage 

Download the software release to have the testing sequences
```
mkdir 3plex;
wget -O v0.1.1-beta.zip https://github.com/molinerisLab/3plex/archive/refs/tags/v0.1.1-beta.zip;
unzip 0.1.1-beta.zip;
```

Then run the test using docker pulled from docher hub.
```
docker run -u `id -u` -it --rm -v $PWD:$PWD imolineris/3plex:v0.0.2 $PWD/test/ssRNA.fa $PWD/test/dsDNA.fa $PWD/test_out/
```

Check the outuput files

```
ls test_out/*/
```
