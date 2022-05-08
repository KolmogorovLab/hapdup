FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, fenderglass@gmail.com

# update and install dependencies
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
    apt-get -y install cmake git make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev tabix libbz2-dev && \
	apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev liblzma-dev libhdf5-dev && \
	apt-get -y install python3-pip python3-virtualenv virtualenv python3-dev && \
	apt-get clean && \
	apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN ln -s /usr/bin/python3 /usr/bin/python
#RUN which python

RUN python3 --version && \
	python3 -m pip install --upgrade pip && \
	python3 -m pip install cython wheel pysam numpy biopython && \
	python3 -m pip uninstall -y enum34 && \
#disables Cuda, but saved ~2Gb of image size
	python3 -m pip install torch==1.9.0+cpu -f https://download.pytorch.org/whl/cpu/torch_stable.html


#build and Flye
WORKDIR /opt
COPY ./submodules/Flye /opt/flye
RUN cd /opt/flye && python3 setup.py install && \
	rm -rf /opt/flye

#RUN cmake --version

# get PEPPER
WORKDIR /opt
COPY ./submodules/pepper /opt/pepper/
#COPY ./submodules/pepper-private /opt/pepper/
RUN cd /opt/pepper && \
    python3 -m pip install . && \
	rm -rf /opt/pepper

#RUN python3 -m pip install --upgrade pip
#RUN rm -rf /opt/pepper/

# install Margin
WORKDIR /opt
COPY ./submodules/margin/ /opt/margin_dir/
RUN cd margin_dir/ && \
	mkdir build && \
	cd build && \
	cmake .. && \
	make -j 20 && \
	cp margin /usr/local/bin && \
	rm -rf /opt/margin_dir

#install the pipeline
WORKDIR /opt
COPY . /opt/hapdup
RUN cd /opt/hapdup && python3 -m pip install . && \
	rm -rf /opt/hapdup

# setup models/configurations
COPY ./pepper_models /opt/pepper_models/
COPY ./submodules/margin/params /opt/margin_params/

ENV PEPPER_MODEL_DIR "/opt/pepper_models"
ENV MARGIN_CONFIG_DIR "/opt/margin_params/phase"

ENV PYTHONUNBUFFERED "1"
