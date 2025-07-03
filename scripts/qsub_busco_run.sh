#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q mem_512_12h,short-sl7,long-sl7
#$ -l virtual_free=10G,h_rt=6:00:00
#$ -o tmp/
#$ -e tmp/

if [ -z $6 ] ; then
	echo "missing arguments:"
	echo "1. input fasta"
	echo "2. run code (e.g. species code)"
	echo "3. code of BUSCO dataset (will download if not already present)."
	echo "4. BUSCO mode (prot for proteins, tran for transcriptome, geno for genomes (UNTESTED)."
	echo "5. output folder (a folder named busco_<run code>_<dataset code> will be created within this folder"
	echo "6. num threads"
	exit
fi

fas=$(readlink -f $1)
spi=$2
ref=$3
mod=$4
out=$5
nth=$6

timestamp=$(date +%s)

# create output folder if needed
mkdir -p ${out}
cd ${out}

# run busco
busco -i ${fas} -o busco_${spi}_${ref} -m ${mod} -l ${ref} --download_path ~/data/busco_databases -c ${nth} -f

