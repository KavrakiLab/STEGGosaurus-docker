FROM ubuntu:noble
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y install \
        cmake \
        cython3 \
        doxygen \
        freeglut3-dev \
        gfortran \
        libeigen3-dev \
        libfftw3-dev \
        libfreetype6-dev \
        libglew-dev \
        libglm-dev \
        libpng-dev \
        libxml2-dev \
        python3-pip \
        swig \
        wget \
        vim

RUN ["/bin/bash", "-c", "wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"]
RUN chmod 0755 Miniforge3-Linux-x86_64.sh
RUN bash Miniforge3-Linux-x86_64.sh -b -p /mamba
ENV PATH "/mamba/bin:$PATH"

COPY ./Docker_image/environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml
RUN echo "source activate STEGG_all" >> ~/.bashrc
ENV PATH "/mamba/envs/STEGG_all/bin:$PATH"

# These are okay outside environment.yml because they install to PATH
# Perhaps better if inside yml?
RUN mamba install -c bioconda smina
RUN mamba install -c bioconda autodock-vina

# This for some reason creates issues with smina, I am not sure why
RUN pip install mpire
RUN pip install DockQ
# RUN pip install pdb-tools

# RUN apt update
# RUN ["/bin/bash", "-c", "wget http://security.ubuntu.com/ubuntu/pool/universe/p/python2.7/python2.7_2.7.18-13ubuntu1.5_amd64.deb"]  
# RUN ["/bin/bash", "-c", "wget http://security.ubuntu.com/ubuntu/pool/universe/p/python2.7/libpython2.7-stdlib_2.7.18-13ubuntu1.5_amd64.deb"]  
# RUN ["/bin/bash", "-c", "wget http://security.ubuntu.com/ubuntu/pool/universe/p/python2.7/python2.7-minimal_2.7.18-13ubuntu1.5_amd64.deb"]  
# RUN ["/bin/bash", "-c", "wget http://security.ubuntu.com/ubuntu/pool/universe/p/python2.7/libpython2.7-minimal_2.7.18-13ubuntu1.5_amd64.deb"]  
# RUN apt install -y \
#     ./libpython2.7-minimal_2.7.18-13ubuntu1.5_amd64.deb \
#     ./libpython2.7-stdlib_2.7.18-13ubuntu1.5_amd64.deb \
#     ./python2.7-minimal_2.7.18-13ubuntu1.5_amd64.deb \
#     ./python2.7_2.7.18-13ubuntu1.5_amd64.deb
# RUN apt install -y python2.7

# RUN conda install -y anaconda::mkl

COPY Ape-Gen2.0-main /home/Ape-Gen2.0-main
COPY haddock3_TCR /home/haddock3_TCR
COPY STEGG_controler /home/STEGG_controler
COPY tfold /home/tfold
COPY T-RECS /home/T-RECS

ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/mamba/envs/STEGG_all/lib/"

# Change dir
WORKDIR /home

ENTRYPOINT ["bash"]

# # Build this Dockerfile into an image by getting to the Dockerfile dir and executing:
# # docker build . -t kavrakilab/stegg
