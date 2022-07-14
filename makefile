VERSION?=v0.0.2
all: docker_build

#triplexator_docker_context/min_len_from_10_to_6.patch:
#	cd triplexator; \
#	git diff 4505bba7b3dc8cf4922d71446c755e91c448673c 9abaef685bd6a4183ebfe689971b87681f803df5 > ../min_len_from_10_to_6.patch
#
#triplexator_patch: min_len_from_10_to_6.patch 
#	cp $< triplexator/
#	cd triplexator; \
#	git apply $<


#triplexator.v1.3.2_l6.tar.gz:
#	docker save triplexator:v1.3.2_l6 | gzip > $@
#
#triplexator_v1.3.2_l6-%.sif:
#	docker run -v /var/run/docker.sock:/var/run/docker.sock -v $${PWD}/singularity_images:/output --privileged -t --rm quay.io/singularity/docker2singularity triplexator:v1.3.2_l6

docker_context/v1.3.3.zip: 
	wget -O $@ https://github.com/Gurado/triplexator/archive/refs/tags/v1.3.3.zip

build: docker_context/v1.3.3.zip docker_context/min_len_from_10_to_6.patch
	docker build -t 3plex:$(VERSION) docker_context

docker_context/conda_environment.yaml:
	#conda activate 3plex_v0.1;\
	#conda env export --from-history > $@
	#this modality do not include version if not explicated during install 
