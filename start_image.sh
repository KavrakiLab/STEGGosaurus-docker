#!/bin/bash
docker run --name benevolent_stegosaurus \
  --gpus all -it \
  --ipc=host --shm-size=1g \
  -e OMPI_MCA_opal_cuda_support=1 \
  kavrakilab/stegg:3.0