FROM ubuntu:mantic
MAINTAINER Ivan Molineris <ivan.molineris@unito.it>
ENV LANG=C LC_ALL=C
RUN mkdir /3plex
WORKDIR /3plex
RUN apt-get update && apt-get install -y curl cmake gcc build-essential libgomp1 libboost-all-dev vim ssh-client wget unzip python3 python3-yaml git gawk && apt-get clean
RUN git clone https://github.com/UDC-GAC/pato
RUN cd pato && make gnu -
RUN mkdir /opt/bin/
RUN mv /3plex/pato/target/gnu/pato.release /opt/bin/pato
COPY 3plex_env.yaml /3plex
SHELL ["/bin/bash", "-c","-l"]
RUN mkdir mamba && export MAMBA_ROOT_PREFIX=/3plex/mamba && cd / && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba && eval "$(./bin/micromamba shell hook -s posix)" && micromamba shell init --shell bash && ./bin/micromamba shell init -s bash -p ~/micromamba
RUN export MAMBA_ROOT_PREFIX=/3plex/mamba && eval "$(micromamba shell hook --shell bash)" && source ~/.bashrc && micromamba activate && micromamba env create -f /3plex/3plex_env.yaml 
COPY . /3plex
RUN mkdir /.cache && chmod -R 777 /3plex && chmod -R 777 /.cache
ENTRYPOINT ["/3plex/3plex.py"]

