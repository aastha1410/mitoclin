#!/bin/bash

echo "Creating environment..."
conda create -n mitoclin python=3.10 -y
conda activate mitoclin

echo "Installing tools..."
conda install -c bioconda bwa samtools bcftools fastqc trim-galore gatk4 haplocheck -y

echo "Setup complete!"
