FROM ubuntu:mantic
MAINTAINER Ivan Molineris <ivan.molineris@unito.it>
ENV LANG=C LC_ALL=C
#TODO: do not copy the actual environment in local/env/3plex
COPY . /3plex
WORKDIR /3plex
RUN apt-get update && apt-get install -y cmake gcc build-essential libgomp1 libboost-all-dev vim ssh-client wget unzip python3 git && apt-get clean
RUN git clone https://github.com/UDC-GAC/pato
RUN cd pato && make gnu -
RUN mkdir /opt/bin/
RUN mv /3plex/pato/target/gnu/pato.release /opt/bin/pato
#RUN wget -O miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#RUN chmod ugo+x miniconda.sh 
#RUN ./miniconda.sh -b -p /opt/conda && \
#    rm miniconda.sh && \
#    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
#    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/profile
#SHELL ["/bin/bash", "-c","-l"]
#RUN conda env create -f 3plex.yaml
ENTRYPOINT ["3plex.py"]


