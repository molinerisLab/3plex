3PLEX_DOCKER_VERSION=v0.0.1
TESTS=ssRNA-dsDNA.tpx ssRNA-dsDNA.tpx.stability raw.tpx.stability

all: $(TESTS)

clean:
	rm $(TESTS)

ssRNA-dsDNA.tpx: ssRNA.fa dsDNA.fa
	docker run -u `id -u` --rm -v $(PWD):$(PWD) 3plex:$(3PLEX_DOCKER_VERSION) bash -c "cd $(PWD);\
		triplexator -l 10  -L -1  -e 10  -g 10  -fr on  -c 3  -fm 0  -of 1   -rm 2   -p 4  \
			 -ss $< \
			 -ds $(word 2,$^) \
			 -o $@;\
			 echo triplexator retstatus: $$?\
	"
	@cmp $@ $@.cfr || echo "TEST FAILD" 

%.tpx.stability: %.tpx
	../docker_context/tpx_thermal_stability.py < $< | sort > $@
	@cmp $@ $@.cfr || echo "TEST FAILD" 

Z85996.1_ss20-dsDNA.tpx.stability.tts_genome_coords.bed: Z85996.1.fa dsDNA.fa dsDNA.bed config.yaml
	docker run --rm -v $(PWD):$(PWD) 3plex:$(3PLEX_DOCKER_VERSION) bash -c "cd $(PWD);\
		source /etc/profile;\
		conda activate 3plex_v0.1;\
		snakemake -j 1 -p $@\
	"
	@cmp $@ $@.cfr || echo "TEST FAILD" 

#TODO need the configuration 
Z85996.1_ss20.fa: Z85996.1.fa
	export PATH=$$PATH:../docker_context/;\
	conda activate 3plex_v0.1;\
	snakemake -s Snakefile_RNA_secondary_structure -j 1 -f Z85996.1_ss20.fa
	@cmp $@ $@.cfr || echo "TEST FAILD" 
