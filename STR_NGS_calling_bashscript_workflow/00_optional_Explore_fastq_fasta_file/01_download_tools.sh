#!/bin/bash

#Install condacolab for conda
pip install -q condacolab
import condacolab
condacolab.install()

# Check the installation
conda --version

# Update conda
conda update conda

#Install sra toolkit from conda
conda install bioconda::sra-tools
