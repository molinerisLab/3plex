TESTS=ssRNA-dsDNA.tpx ssRNA-dsDNA.tpx.stability raw.tpx.stability

all: $(TESTS)

clean:
	rm $(TESTS)

ssRNA-dsDNA.tpx: ssRNA.fa dsDNA.fa
	docker run -u `id -u` --rm -v $(PWD):$(PWD) triplexator:v1.3.2_l6 bash -c "cd $(PWD);\
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

environment.yaml:
	conda activate 3plex_v0.1;\
	conda env export --from-history > $@
