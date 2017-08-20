#!/bin/bash
#Create Snakemake Environment
#
#installed Miniconda 3 http://snakemake.readthedocs.io/en/stable/tutorial/setup.html
# wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh
conda env create --name loqsproject --file environment.yaml
