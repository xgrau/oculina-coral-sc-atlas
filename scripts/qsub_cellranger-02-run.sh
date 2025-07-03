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
cri=$(readlink --canonicalize $2) # path to cellranger index
fqf=$(readlink --canonicalize $3) # path to folder with cell ranger-formatted fastq files
oum=$(readlink --canonicalize $4) # output folder for mapping data
nth=$5

# path to cellranger
cellranger_path="/home/xavi/Programes/cellranger-7.2.0/cellranger"

# mapping, forcing to 10k cells
fqf=$(readlink --canonicalize ${fqf})
rid=$(basename ${fqf})
echo $fqf
echo "cellranger | mapping cells, species=${sid} sample=${rid}..."
cd ${oum}
${cellranger_path} \
	count \
	--id=map_${sid}_${rid}_xkc \
	--sample=${rid} \
	--fastqs=${fqf} \
	--transcriptome=${cri} \
	--jobmode 'local' \
	--localmem=120 \
	--chemistry 'auto' \
	--localcores ${nth}
cd ${oui}
