#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q short-sl7,mem_512_12h,long-sl7
#$ -l virtual_free=20G,h_rt=6:00:00
#$ -o tmp/
#$ -e tmp/

# input
sid=$1 # sps id
fas=$(readlink --canonicalize $2) # genome fasta file (can be gzipped)
gtf=$(readlink --canonicalize $3) # GTF
fql=$(readlink --canonicalize $4) # file with a list of paths to cell ranger-formatted fastq files
our=$(readlink --canonicalize $5) # output folder for reference genome
oum=$(readlink --canonicalize $6) # output folder for mapping data
nth=$7

# retain
timestamp=$(date +%s)
oui=$(pwd)

# path to cellranger
cellranger_path="/home/xavi/Programes/cellranger-6.1.1/cellranger"

# create dirs
mkdir -p ${our}
mkdir -p ${oum}

# create cellranger index
# this job will run first, and the mapping jobs will remain on hold until this is ready
echo "cellranger | mkref ${sid} start"
bash ../scripts/qsub_cellranger-01-db.sh \
		${sid} \
		${fas} \
		${gtf} \
		${our} \
		${nth}

# run cellranger mapping steps
# run one job per read pair. Jobs will remain on hold until index is ready.
echo "cellranger | mapping 10k cells, ${sid} launch jobs and put on hold"
i=0
while read fqf ; do
	let i=i+1
	bash ../scripts/qsub_cellranger-02-run.sh \
			${sid} \
			${our}/gdb_${sid} \
			${fqf} \
			${oum} \
			${nth}
done < ${fql}

# finish
echo "cellranger | ${i} mapping jobs submitted for ${sid}"
