all: docker_build

triplexator-1.3.2-Linux.tar.gz:
	wget -O $@ https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/triplexator/triplexator-1.3.2-Linux.tar.gz

triplexator-1.3.2_l6-Linux.tar.gz: triplexator
	tar -cvzf $@ $<

triplexator_docker_context/min_len_from_10_to_6.patch:
	cd triplexator; \
	git diff 4505bba7b3dc8cf4922d71446c755e91c448673c 9abaef685bd6a4183ebfe689971b87681f803df5 > ../min_len_from_10_to_6.patch

#triplexator_patch: min_len_from_10_to_6.patch 
#	cp $< triplexator/
#	cd triplexator; \
#	git apply $<


#triplexator_docker_context/%.tar.gz: %.tar.gz
#	ln  $< $@


triplexator_docker_context/v1.3.3.zip: 
	wget -O $@ https://github.com/Gurado/triplexator/archive/refs/tags/v1.3.3.zip

docker_build: triplexator_docker_context/v1.3.3.zip triplexator_docker_context/min_len_from_10_to_6.patch
	docker build -t triplexator:v1.3.2_l6 triplexator_docker_context

triplexator.v1.3.2_l6.tar.gz:
	docker save triplexator:v1.3.2_l6 | gzip > $@

triplexator_v1.3.2_l6-%.sif:
	docker run -v /var/run/docker.sock:/var/run/docker.sock -v $${PWD}/singularity_images:/output --privileged -t --rm quay.io/singularity/docker2singularity triplexator:v1.3.2_l6

