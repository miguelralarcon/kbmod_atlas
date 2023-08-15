FROM nvidia/cuda:12.0.0-base-ubuntu20.04
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for building
RUN apt-get update && apt-get install -y cmake \
    python3.8 python3.8-dev python3-pip wget git \
    build-essential cuda-toolkit-12-0\
    && rm -rf /var/lib/apt/lists/*
RUN python3.8 -m pip install --upgrade pip
RUN python3.8 -m pip install pipenv
WORKDIR /app

# Copy the code from the local directory to the container
ADD . .

# Build the project
RUN python3.8 -m pip install .
RUN python3.8 -m pip install .[analysis]
RUN python3.8 -m pip install cupy-cuda12x

# Install JupyterLab
RUN pip3 install jupyterlab

# Give execution permissions to warpsky
RUN chmod +x parsers/warpsky

# Set the entry point command to launch JupyterLab
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--allow-root"]


# Docker build
# docker build -t kbmod .

# Docker run
# sudo docker run --restart=always -d -p 8187:8888 --runtime nvidia -e JUPYTER_TOKEN=dirac -e JUPYTER_ENABLE_LAB=yes -v /mnt/vast/atlas-tdo:/app/images --name kbmod kbmod
