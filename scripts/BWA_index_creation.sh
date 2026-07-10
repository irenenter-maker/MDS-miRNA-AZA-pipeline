#/bin/bash

mkdir -p data/mature_hsa_index
mkdir -p data/hairpin_hsa_index

bwa index -a is -p data/mature_hsa_index/mature_hsa data/mature_hsa.fa
bwa index -a is -p data/hairpin_hsa_index/hairpin_hsa data/hairpin_hsa.fa