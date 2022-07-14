FROM debian/eol:wheezy
MAINTAINER Ivan Molineris <ivan.molineris@unito.it>

ENV LANG=C LC_ALL=C
RUN apt-get update && apt-get install -y cmake g++ build-essential libgomp1 libboost-all-dev vim ssh-client wget unzip python3 git && apt-get clean
WORKDIR /tmp
#RUN git clone git@github.com:Gurado/triplexator.git
ADD v1.3.3.zip /tmp/
RUN unzip v1.3.3.zip
WORKDIR /tmp/triplexator-1.3.3
ADD min_len_from_10_to_6.patch /tmp/triplexator-1.3.3
RUN git apply min_len_from_10_to_6.patch 
RUN mkdir -p build/Release && cd build/Release && cmake ../.. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -G "Unix Makefiles" && make
RUN cp bin/triplexator /usr/bin/


ENV MINICONDA_VERSION="4.7.12.1"
RUN wget -O miniconda.sh https://repo.anaconda.com/miniconda/Miniconda2-${MINICONDA_VERSION}-Linux-x86_64.sh 
RUN chmod ugo+x miniconda.sh 
RUN ./miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/profile
SHELL ["/bin/bash", "-c","-l"]
WORKDIR /tmp/
ENV CONDA_ENV_VERSION="0.2"
ADD conda_environment.yaml /tmp/conda_environment.yaml
RUN conda env create -f conda_environment.yaml

ADD distMedian.py /usr/bin/distMedian.py
ADD tpx_thermal_stability.py /usr/bin/tpx_thermal_stability.py
ADD find_best /usr/bin/find_best
RUN mkdir -p /opt/local/lib/perl
ADD MakefileParser.pm /opt/local/lib/perl
RUN apt-get install -y gawk
ADD bawk /usr/bin/bawk
ADD fasta2oneline /usr/bin/fasta2oneline
ADD bsort /usr/bin/bsort
ADD translate /usr/bin/translate
ADD id2count /usr/bin/id2count
#move upstairs 
