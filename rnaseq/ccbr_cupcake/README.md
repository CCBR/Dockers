# cDNA_Cupcake

Creating a virtual environment for cDNA_Cupcake.

### One-time site installation

The steps directly below only need to be run once to install Cupcake on your target system.
Once installed the virtual enviroment will only need to be activated prior to running
any of Cupcake’s commands.

```bash
# Cupcake requires python >= 3.7
module load python/3.8

# Initialize virtual environment
python -m venv .env

# Activate the env
source .env/bin/activate

# Upgrade the to the latest version of pip
pip install --upgrade pip

# Clone the github repo
git clone https://github.com/Magdoll/cDNA_Cupcake.git
cd cDNA_Cupcake/

# Install requirement packages for setup.py
pip install -r requirements.txt
pip install Cython

# Install cDNA_Cupcake
python setup.py build
python setup.py install
PATH=$PATH:${PWD}/sequence/
```

### Run cupcake after installing on target system

As mentioned above, the steps above only need to be run once. This is to install
cupcake on your target system. Once cupcake is install, a user can activate and run
cupcake using the following steps:

```bash
# If the virtual env is not activated
# assumes you are in the same directory as the virtual environment
# If not in the same dir, prefix the activate with the absolute PATH to the virtual
# environment. (i.e. source /path/to/virtual/environment/.env/bin/activate)
source .env/bin/activate

# You may need to add this to your $PATH: sequence/ directory in cDNA_Cupcake repo
export PATH=$PATH:/path/to/github/repo/install/cDNA_Cupcake/sequence

# Run whatever cupcake commands you need for your analysis
<run your cupcake commands>

deactivate
```
### Build from Dockerfile

The instructions above are for site installations of cDNA_cupcake. The instructions below are for building a Docker image from the provided Dockerfile. 

```bash
# See listing of images on computer
docker image ls

# Build
docker build --tag=ccbr_cupcake:v0.0.1 .

# Updating tag(s) before pushing to DockerHub
docker tag ccbr_cupcake:v0.0.1 skchronicles/ccbr_cupcake:v0.0.1
docker tag ccbr_cupcake:v0.0.1 skchronicles/ccbr_cupcake        # latest
docker tag ccbr_cupcake:v0.0.1 nciccbr/ccbr_cupcake:v0.0.1
docker tag ccbr_cupcake:v0.0.1 nciccbr/ccbr_cupcake             # latest

# Check out new tag(s)
docker image ls

# Peak around the container: verify things run correctly
docker run -ti ccbr_cupcake:v0.0.1 /bin/bash

# Push new tagged image to DockerHub
docker push skchronicles/ccbr_cupcake:v0.0.1
docker push skchronicles/ccbr_cupcake:latest
docker push nciccbr/ccbr_cupcake:v0.0.1
docker push nciccbr/ccbr_cupcake:latest
```
