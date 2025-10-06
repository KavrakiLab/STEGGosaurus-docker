#!/bin/bash

# # If CUDA is needed, install cuda toolkit in the container
# mamba activate STEGG_all
# mamba install -y -c conda-forge \
#   cuda-toolkit=12.8 \
#   cuda-nvcc=12.8 \
#   cuda-cudart-dev=12.8
# export CUDA_HOME="$CONDA_PREFIX"
# export PATH="$CUDA_HOME/bin:$PATH"
# export LD_LIBRARY_PATH="$CUDA_HOME/lib:$CUDA_HOME/lib64:${LD_LIBRARY_PATH}"

# Run the docker container with GPU support
# docker run --name benevolent_stegosaurus_gpu \
#   --gpus all -it \
#   --ipc=host --shm-size=1g \
#   -e OMPI_MCA_opal_cuda_support=1 \
#   kavrakilab/stegg:latest

# Run the docker container without GPU support
docker run --name benevolent_stegosaurus_no_gpu -it --ipc=host --shm-size=1g kavrakilab/stegg:latest