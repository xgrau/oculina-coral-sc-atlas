#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7
#$ -l virtual_free=50G,h_rt=720:00:00
#$ -o tmp/
#$ -e tmp/

# input
fas=$1 # folder with fasta input
nth=$2 # num threads

python ~/Programes/Broccoli/broccoli.py -threads ${nth} -dir ${fas} -nb_hits 10 -phylogenies ml -kmer_size 10000
