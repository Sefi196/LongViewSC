#!/bin/bash

# Set the environment name
ENV_NAME="LongViewSC_env"

# Step 1: Create the Conda environment from the provided yml file
echo "Creating Conda environment from environment.yml"
conda env create -f environment.yml

# Step 2: Activate the Conda environment
echo "Activating the Conda environment: $ENV_NAME"
conda activate $ENV_NAME

# Step 3: Install R dependencies using devtools
echo "Installing ggtranscript R package..."
R -e "devtools::install_github('dzhang32/ggtranscript')"

# Step 4: Confirm successful installation
echo "Installation complete! Your environment is set up and the ggtranscript package has been installed."
