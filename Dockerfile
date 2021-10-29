FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, fenderglass@gmail.com

# update and install dependencies
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
    apt-get -y install wget git cmake make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev tabix libbz2-dev && \
	apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev liblzma-dev libhdf5-dev && \
	apt-get -y install python3-pip python3-virtualenv virtualenv && \
	apt-get clean && \
	apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && \
	apt-get -y install time git make wget autoconf gcc g++ && \
	apt-get -y install autoconf bzip2 lzma-dev zlib1g-dev && \
	apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev && \
	apt-get -y install liblzma-dev libhdf5-dev libncurses5-dev libncursesw5-dev && \
	apt-get -y install python3-dev python3-pip && \
	apt-get -y install cmake && \
	apt-get -y install protobuf-compiler libprotobuf-dev && \
	apt-get clean && \
	apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN ln -s /usr/bin/python3 /usr/bin/python
RUN which python

RUN python3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install cython
RUN python3 -m pip install wheel
RUN python3 -m pip uninstall -y enum34

RUN python3 -m pip install pysam
RUN python3 -m pip install numpy

#build and Flye
WORKDIR /opt
COPY ./submodules/Flye /opt/flye
RUN cd /opt/flye && python3 setup.py install
RUN rm -rf /opt/flye

RUN cmake --version

# get PEPPER
WORKDIR /opt
#COPY ./submodules/pepper /opt/pepper/
COPY ./submodules/pepper-private /opt/pepper/
RUN cd /opt/pepper && \
    python3 -m pip install .

RUN python3 -m pip install --upgrade pip
RUN rm -rf /opt/pepper/

# install Margin
WORKDIR /opt
COPY ./submodules/margin/ /opt/margin_dir/
RUN cd margin_dir/ && \
	mkdir build && \
	cd build && \
	cmake .. && \
	make -j 20 && \
	cp margin /usr/local/bin
RUN rm -rf /opt/margin_dir

#install the pipeline
WORKDIR /opt
COPY . /opt/hapdup
RUN cd /opt/hapdup && python3 -m pip install .
RUN rm -rf /opt/hapdup

# setup models/configurations
COPY ./pepper_models /opt/pepper_models/
COPY ./submodules/margin/params /opt/margin_params/

ENV PEPPER_MODEL "/opt/pepper_models/PEPPER_VARIANT_ONT_R941_GUPPY5_SUP_V6.pkl"
ENV MARGIN_MODEL "/opt/margin_params/phase/allParams.haplotag.ont-r94g507.incSupAln.json"

ENV PYTHONUNBUFFERED "1"
