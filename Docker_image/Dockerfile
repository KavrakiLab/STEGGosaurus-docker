# syntax=docker/dockerfile:1.6

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
        vim && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    chmod 0755 miniforge.sh && \
    bash miniforge.sh -b -p /mamba && \
    rm miniforge.sh
ENV PATH="/mamba/bin:$PATH"

COPY ./Docker_image/environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy
# RUN echo "source activate STEGG_all" >> ~/.bashrc
SHELL ["/bin/bash", "-lc"]
RUN mamba shell init --shell bash
RUN printf '%s\n' \
  'eval "$(mamba shell hook --shell bash)"' \
  'mamba activate STEGG_all' \
  > /etc/profile.d/activate_stegg.sh
RUN eval "$(mamba shell hook --shell bash)" \
 && mamba activate STEGG_all \
 && python -c "import sys; print(sys.executable)"

ENV PATH="/mamba/envs/STEGG_all/bin:$PATH"

# These are okay outside environment.yml because they install to PATH
# Perhaps better if inside yml?
RUN mamba install -c bioconda smina && mamba clean -afy
RUN mamba install -c bioconda autodock-vina && mamba clean -afy

# This for some reason creates issues with smina, I am not sure why
RUN pip install mpire DockQ && pip cache purge
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

# Save pytorch model files
ENV TORCH_HOME=/opt/torch_models
RUN --mount=type=bind,source=models,target=/__models,readonly \
    mkdir -p "$TORCH_HOME/hub/checkpoints" && \
    cp -a /__models/. "$TORCH_HOME/hub/checkpoints/" && \
    chmod -R 755 "$TORCH_HOME"

# Copy required files
RUN --mount=type=bind,source=.,target=/__ctx,readonly \
    cp -a \
      /__ctx/Ape-Gen2.0-main \
      /__ctx/haddock3_TCR \
      /__ctx/STEGG_controler \
      /__ctx/tfold \
      /__ctx/T-RECS \
      /__ctx/RCD_required_files \
      /__ctx/FASPR_required_files \
      /home/ && \
    chmod -R 755 /home

ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/mamba/envs/STEGG_all/lib/"

# Change dir
WORKDIR /home/STEGG_controler

ENTRYPOINT ["bash"]

# # Build this Dockerfile into an image by getting to the Dockerfile dir and executing:
# # docker build . -t kavrakilab/stegg
