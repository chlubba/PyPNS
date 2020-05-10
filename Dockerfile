# Distributed under the terms of the Modified BSD License.
ARG BASE_CONTAINER=jupyter/scipy-notebook
FROM $BASE_CONTAINER
LABEL maintainer="Russell <rjjarvis@asu.edu>"
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    fonts-dejavu \
    gfortran \
    gcc && \
    rm -rf /var/lib/apt/lists/*
RUN apt-get upgrade    
RUN apt-get install -y ncurses-base ncurses-bin
RUN apt-get update && apt-get install -y bzip2 ca-certificates automake libtool  \
                       libncurses5-dev libreadline-dev libgsl0-dev cmake ssh
USER jovyan

WORKDIR $HOME
RUN \
  wget http://www.neuron.yale.edu/ftp/neuron/versions/v7.7/nrn-7.7.tar.gz && \
  tar -xzf nrn-7.7.tar.gz && \
  rm nrn-7.7.tar.gz
WORKDIR $HOME/nrn-7.7
ENV PATH /usr/bin/python3/python:/opt/conda/bin:/opt/conda/bin/conda:/opt/conda/bin/python:$PATH
RUN ./configure --prefix=`pwd` --without-iv --with-nrnpython=/opt/conda/bin/python3
USER root
RUN sudo make all && \
     make install
USER jovyan
WORKDIR src/nrnpython
RUN python setup.py install
RUN python -c "import neuron"
ENV NEURON_HOME $HOME/nrn-7.7/x86_64
ENV PATH $NEURON_HOME/bin:$PATH
WORKDIR $HOME/work/extra_work
WORKDIR $HOME/work
RUN pip install --upgrade matplotlib
# RUN conda install tk
# RUN pip install tkinter
RUN python -c "import tkinter"
#import matplotlib as mpl;mpl.use('TkAgg')"
RUN python -c "import tkinter"
#; import tk"
RUN git clone https://github.com/fun-zoological-computing/PyPNS
WORKDIR PyPNS
RUN pip install -e .

WORKDIR mods
RUN nrnivmodl
#RUN pip install tk
#RUN pip install --upgrade tk


RUN python -c "import neuron"
RUN python -c "import PyPNS"
#RUN python -c "import matplotlib as mpl;mpl.use('TkAgg'); import PyPNS"
WORKDIR $HOME/work/PyPNS/
RUN ls mods/*
RUN cp mods/*.mod .
RUN nrnivmodl
RUN python test.py
WORKDIR $HOME
#ENTRYPOINT /bin/bash