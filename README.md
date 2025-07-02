# Comparative single-cell omic analyses of the stony coral *Oculina patagonica* <!-- omit from toc -->

- [Introduction](#introduction)
- [Accessing the data](#accessing-the-data)
- [Reproducing the analyses](#reproducing-the-analyses)
	- [Genome annotation](#genome-annotation)
		- [*de novo* gene annotation of *Oculina patagonica*](#de-novo-gene-annotation-of-oculina-patagonica)
		- [K-mer frequency plots of *Oculina* assemblies](#k-mer-frequency-plots-of-oculina-assemblies)
	- [Comparative genomic analyses](#comparative-genomic-analyses)
		- [Orthology analyses](#orthology-analyses)
		- [Gene family annotation](#gene-family-annotation)
		- [Gene ontology mappings](#gene-ontology-mappings)
		- [Species tree of anthozoans](#species-tree-of-anthozoans)
		- [Evolution of gene family gain/loss/expansion](#evolution-of-gene-family-gainlossexpansion)
		- [Analysis of tandem duplications](#analysis-of-tandem-duplications)
		- [Microsynteny analyses](#microsynteny-analyses)
		- [Macrosynteny analyses](#macrosynteny-analyses)
		- [Whole-genome alignment and conservation](#whole-genome-alignment-and-conservation)
		- [Ks rate analysis of whole-genome duplications](#ks-rate-analysis-of-whole-genome-duplications)
		- [Transposon complement analysis](#transposon-complement-analysis)
	- [Single-cell transcriptomic atlas analyses](#single-cell-transcriptomic-atlas-analyses)
		- [Mapping scRNAseq to corals](#mapping-scrnaseq-to-corals)
		- [Cluster and annotate scRNAseq atlases](#cluster-and-annotate-scrnaseq-atlases)
	- [Cross-species single-cell transcriptome analyses](#cross-species-single-cell-transcriptome-analyses)
		- [Cell type trees across species](#cell-type-trees-across-species)
		- [Use SAMap to align cell types across species](#use-samap-to-align-cell-types-across-species)

## Introduction

Folder organisation:

* `data/`: reference genomes, sequencing data, scRNA-seq data.
* `results_scatlas/`: scRNA-seq atlases.
* `results_annotation/`: genome annotation and comparative analyses.

Each species is referred to with a short acronym all through this project. For the three scleractinian coral single-cell atlases here produced, we use:

* `Ocupat` for *Oculina patagonica*
* `Amil` for *Acropora millepora*
* `Spin` for *Stylophora pistillata*

We have compared these three corals to other species, namely:

* `Ocuarb` for *Oculina arbuscula*, another scleractinian coral
* `Xesp` for *Xenia* sp., an octocorallian soft coral
* `Nvec` for *Nematostella vectensis*, a sea anemone

The complete list of species used in the comparative transcriptomic and genomic analyses, with the relevant data sources, can be found [here](data/taxon_sampling.md).

## Accessing the data

We provide the following [raw data](results_scatlas/mapping/scdb_seurat) for each coral species:

* Seurat objects
* Metacell annotation tables
* Cell type annotation tables
* Footprints

If you have any queries, feel free to let me know in the [Issues section](https://github.com/xgrau/oculina-coral-sc-atlas/issues).

## Reproducing the analyses

### Genome annotation

#### *de novo* gene annotation of *Oculina patagonica*

```bash
## TODO ##
```

#### K-mer frequency plots of *Oculina* assemblies

We use [`KAT`](https://kat.readthedocs.io/en/latest/walkthrough.html#heterozygous-genomes) to create k-mer frequency plots for the raw and processed assemblies.

To run from the `data` folder.

```bash
# cdir
cd sequencing_nanopore_Oculina

# processed assembly
kat comp -t 12 concatenated_nanopore_reads.fastq.gz ../reference/Ocupat_gDNA.fasta -o kat_comp_raw --output_type=pdf --output_hists -v
kat plot spectra-cn kat_comp-main.mx --x_max 100 --output_type=pdf -o kat_comp-main.spectra.pdf
kat plot density kat_comp-main.mx --x_max 100 --output_type=pdf -o kat_comp-main.density.pdf
# raw assembly
kat comp -t 36 concatenated_nanopore_reads.fastq.gz intermediate_assemblies_flye01/assembly.fasta -o kat_comp_raw --output_type=pdf --output_hists -v
kat plot spectra-cn kat_comp_raw-main.mx --x_max 100 --output_type=pdf -o kat_comp_raw-main.spectra.pdf
kat plot density kat_comp_raw-main.mx --x_max 100 --output_type=pdf -o kat_comp_raw-main.density.pdf
```

### Comparative genomic analyses

#### Orthology analyses

To launch from the `data/` folder.

1. Anthozoa-level orthology database (used for intra-anthozoan comparative analyses):

```bash
mkdir -p orthology_Anthozoa/dataset
for i in Nvec Scocal Actieq Metsen Dialin Exapal Adig Amil Gfas Fspp Gasp Spis Ocupat Ocuarb Pocdam Xesp Dgig ; do
gffread reference/${i}_long.annot.gtf -g reference/${i}_gDNA.fasta -y orthology_Anthozoa/dataset/${i}_long.pep.fasta
bioawk -c fastx '{ n=gsub(/\./, "X", $2) ; print ">"$1"\n"$2 }' orthology_Anthozoa/dataset/${i}_long.pep.fasta > orthology_Anthozoa/dataset/${i}_long.pep.fasta.tmp && mv orthology_Anthozoa/dataset/${i}_long.pep.fasta.tmp orthology_Anthozoa/dataset/${i}.fasta && rm orthology_Anthozoa/dataset/${i}_long.pep.fasta
done
ls orthology_Anthozoa/dataset
# launch broccoli
cd orthology_Anthozoa/
bash ../../scripts/qsub_broccoli-ml.sh dataset 6

# run busco for each species
conda activate busco
mkdir -p reference/busco
for i in Nvec Scocal Actieq Metsen Dialin Exapal Adig Amil Gfas Fspp Gasp Spis Ocupat Ocuarb Pocdam Xesp Dgig Acrcer Acrpal Porlut Pocver; do
bash ../scripts/qsub_busco_run.sh reference/${i}_long.pep.fasta ${i} metazoa_odb10 prot reference/busco 4
done
```

1. Metazoa-level orthology database (used to obtain named orthologs for selected anthozoans):

```bash
mkdir -p orthology_Metazoa/dataset
for i in Nvec Exapal Amil Spis Ocupat Xesp Dgig Hvul Turdoh Chem Aaur Rhoesc Mmus Hsap Bralan Spur Astrub Skow Dmel Pricau Owefus Lgig Tadh Hhon ; do
gffread ../../genomes/data/${i}_long.annot.gtf -g ../../genomes/data/${i}_gDNA.fasta -y orthology_Metazoa/dataset/${i}_long.pep.fasta
bioawk -c fastx '{ n=gsub(/\./, "X", $2) ; print ">"$1"\n"$2 }' orthology_Metazoa/dataset/${i}_long.pep.fasta > orthology_Metazoa/dataset/${i}_long.pep.fasta.tmp && mv orthology_Metazoa/dataset/${i}_long.pep.fasta.tmp orthology_Metazoa/dataset/${i}.fasta && rm orthology_Metazoa/dataset/${i}_long.pep.fasta
done
ll -h orthology_Metazoa/dataset
head orthology_Metazoa/dataset/Ocupat.fasta
# launch broccoli
cd orthology_Metazoa/
bash ../../scripts/qsub_broccoli-ml.sh dataset 24

# run busco for each species
conda activate busco
mkdir -p reference/busco
for i in Dgig Hvul Turdoh Chem Aaur Rhoesc Mmus Hsap Bralan Spur Astrub Skow Dmel Pricau Owefus Lgig Tadh Hhon; do
bash ../scripts/qsub_busco_run.sh reference/${i}_long.pep.fasta ${i} metazoa_odb10 prot reference/busco 4
done
```

1. Consolidate orthogroups. Beware, this step requires having run dedicated phylogenies to annotate transcription factors (see below: [Gene family annotation](#gene-family-annotation))

```bash
# consolidate OGs and pairs
Rscript s01_prepare_orthology_2023-01-20.R

# ancestral reconstruction with Dollo via Possvm
python ../scripts/possvm_reconstruction.py -tree species_tree.Anthozoa.newick -ort orthology_Anthozoa/orthogroup_conservation.csv -out orthology_Anthozoa/orthogroup_conservation.possvm
python ../scripts/possvm_reconstruction.py -tree species_tree.Metazoa.newick  -ort orthology_Metazoa/orthogroup_conservation.csv  -out orthology_Metazoa/orthogroup_conservation.possvm
# with Ocuarb:
python ../scripts/possvm_reconstruction.py -tree species_tree.Anthozoa_plus.newick -ort orthology_Anthozoa_plus/orthogroup_conservation.csv -out orthology_Anthozoa_plus/orthogroup_conservation.possvm
python ../scripts/possvm_reconstruction.py -tree species_tree.Metazoa_plus.newick -ort orthology_Metazoa_plus/orthogroup_conservation.csv -out orthology_Metazoa_plus/orthogroup_conservation.possvm
```

#### Gene family annotation

To run from the `results_annotation` subfolder.

1. Launch phylogenies

```bash
# concatenate fasta
mkdir -p data
cat ../data/orthology_Metazoa/dataset/*.fasta > data/seq_Metazoa_2024-05-14.fasta
esl-sfetch --index data/seq_Metazoa_2024-05-14.fasta

# First, launch dedicated phylogenies
while read -a i ; do bash s01_hmmsearches_v10_2021-12-01.sh ${i[0]} ${i[1]} ${i[2]} ${i[3]} data/seq_Metazoa_2024-05-14.fasta Y tfs results_alignments/ results_searches/ ${i[4]} ; done < <(grep -w TF data/gene_families_searchinfo.csv)

# Second, parse phylogenies with Possvm
# Get mouse gene name for reference
wget http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz
gunzip Mus_musculus.GRCm38.pep.all.fa.gz
paste <(grep -o "transcript:[^ ]*" Mus_musculus.GRCm38.pep.all.fa | sed "s/transcript:/Mmus_/" | sed "s/\.[0-9]*$//") <(grep -o "gene_symbol:[^ ]*" Mus_musculus.GRCm38.pep.all.fa | sed "s/gene_symbol://") > data/mus_gene_names.csv
rm Mus_musculus.GRCm38.pep.all.fa
# copy gene trees
mkdir -p results_gene_trees
for set in tfs ; do
cp results_alignments/${set}.*.treefile results_gene_trees
done

# Run possvm
for set in tfs ; do
for i in results_gene_trees/${set}.*.treefile ; do 
possvm -i $i -o results_gene_trees/ -r data/mus_gene_names.csv -refsps Mmus -itermidroot 10 -min_support_transfer 50 -cut_gene_names 100 -outgroup data/outgroup_species.txt -ogprefix $(basename $i | cut -f2,3 -d '.'). -p $(basename $i | sed "s/.seqs.iqtree.treefile/.possvm/")
done
done
```

2. Gene lists based on Pfam domain presence (comprehensive and species-specific, they don't depend on gene trees):

```bash
mkdir -p results_gene_lists
for i in Amil Ocupat Spis Xesp Nvec Ocuarb ; do 
pfam_doms=../data/reference/${i}_long.pep.pfamscan_archs.csv
ls ${pfam_doms}
gefamclass=$(cut -f6 data/gene_families_searchinfo.csv | sort -u)
for cla in ${gefamclass} ; do
echo "$i -> $cla"
gefamid=$(fgrep -w ${cla} data/gene_families_searchinfo.csv | cut -f1)
for gef in ${gefamid} ; do
	fgrep -f <(fgrep -w ${gef} data/gene_families_searchinfo.csv | cut -f2 | tr ',' '\n') ${pfam_doms} | cut -f1 | sed "s/$/\t${gef}/"
done | sort -k1,2 | awk '{ if(a[$1]) a[$1]=a[$1]","$2 ; else a[$1]=$2 ; }END { for (i in a) { print i, a[i]; } }' > results_gene_lists/${i}_list_transcripts.${cla}.txt
awk 'NR==FNR { l[$1]=$2;next} { for(i = 1; i <= NF; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' ../data/reference/${i}_all.annot.pep.gene2tx.txt results_gene_lists/${i}_list_transcripts.${cla}.txt > results_gene_lists/${i}_list_genes.${cla}.txt
done
done
```

3. Gene lists for signal peptides, transmembrane domains and, by inference, secreted proteins:

```bash
# transmembrane domains
for i in Amil Ocupat Spis Xesp Nvec Ocuarb; do 
fgrep -w TMhelix ../data/reference/${i}_long.pep.tmhmm.txt | cut -f1 > results_gene_lists/${i}_list_transcripts.with_TM_domain.txt
awk 'NR==FNR { l[$1]=$2;next} { for(i = 1; i <= NF; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' ../data/reference/${i}_all.annot.pep.gene2tx.txt results_gene_lists/${i}_list_transcripts.with_TM_domain.txt > results_gene_lists/${i}_list_genes.with_TM_domain.txt
done
# signal peptides
for i in Amil Ocupat Spis Xesp Nvec Ocuarb; do 
fgrep -w "SP(Sec/SPI)" ../data/reference/${i}_long.pep.signalp.txt | cut -f1 > results_gene_lists/${i}_list_transcripts.with_signalpept.txt
awk 'NR==FNR { l[$1]=$2;next} { for(i = 1; i <= NF; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' ../data/reference/${i}_all.annot.pep.gene2tx.txt results_gene_lists/${i}_list_transcripts.with_signalpept.txt > results_gene_lists/${i}_list_genes.with_signalpept.txt
done
# secreted: signal peptides and no TM domain
for i in Amil Ocupat Spis Xesp Nvec Ocuarb; do 
comm -23 <(sort results_gene_lists/${i}_list_genes.with_signalpept.txt) <(sort results_gene_lists/${i}_list_genes.with_TM_domain.txt) > results_gene_lists/${i}_list_genes.secreted.txt
comm -23 <(sort results_gene_lists/${i}_list_transcripts.with_signalpept.txt) <(sort results_gene_lists/${i}_list_transcripts.with_TM_domain.txt) > results_gene_lists/${i}_list_transcripts.secreted.txt
done
```

4. Gene lists for ribosomal proteins and histones

```bash
# ribosomal
for i in Ocuarb ; do
grep "Ribosomal" ../data/reference/${i}_long.pep.pfamscan_archs.csv | awk 'BEGIN { FS="\t" ; OFS="\t" } { print "ribosomal", $1, $2}' > results_gene_lists/${i}_list_transcripts.ribosomal.txt
awk 'BEGIN { FS="\t" ; OFS="\t" } NR==FNR { l[$1]=$2;next} { for(i = 1; i <= NF; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' ../data/reference/${i}_all.annot.pep.gene2tx.txt results_gene_lists/${i}_list_transcripts.ribosomal.txt > results_gene_lists/${i}_list_genes.ribosomal.txt
done
# histone
for i in Ocuarb ; do
grep "Histone" ../data/reference/${i}_long.pep.pfamscan_archs.csv | awk 'BEGIN { FS="\t" ; OFS="\t" } { print "histones", $1, $2}' > results_gene_lists/${i}_list_transcripts.histones.txt
awk 'BEGIN { FS="\t" ; OFS="\t" } NR==FNR { l[$1]=$2;next} { for(i = 1; i <= NF; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' ../data/reference/${i}_all.annot.pep.gene2tx.txt results_gene_lists/${i}_list_transcripts.histones.txt > results_gene_lists/${i}_list_genes.histones.txt
done
```

#### Gene ontology mappings

From `data` folder:

1. Download Mus data from MIG:

```bash
# download GOs for mous from MIG database and GO database, 2022-11-03 release
wget http://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt -O reference/Mmus_MRK_Sequence_ID_mapping.csv
wget http://current.geneontology.org/annotations/mgi.gaf.gz -O reference/Mmus_MGI.gaf.gz
# MGI to Ensembl mapping
awk 'BEGIN { FS="\t" } { print $1,$13 }' reference/Mmus_MRK_Sequence_ID_mapping.csv | tr '|' ' ' | awk 'NR > 1 { for(i=2; i<=NF ;i++)  { print $1,$i } }' | tr ' ' '\t' | sed "s/\t/\tMmus_/"> reference/Mmus_MRK_Sequence_ID_mapping.dict.txt
fgrep -f <(cut -f1 reference/Mmus_long.pep.pfamscan_archs.csv) reference/Mmus_MRK_Sequence_ID_mapping.dict.txt > reference/Mmus_MRK_Sequence_ID_mapping.dict_long.txt
# GO mappings
gunzip reference/Mmus_MGI.gaf.gz
grep -v "^!" reference/Mmus_MGI.gaf | awk '{ print $2,$5 }' | sort -u | awk '{ if(a[$1]) a[$1]=a[$1]","$2 ; else a[$1]=$2 ; }END { for (i in a) { print i, a[i]; } }' | tr ' ' '\t' > reference/Mmus_MGI.gaf.csv
awk 'NR==FNR { l[$1]=$2;next} { for(i = 1; i <= NF; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' reference/Mmus_MRK_Sequence_ID_mapping.dict_long.txt reference/Mmus_MGI.gaf.csv | grep "^Mmus_" | tr ' ' '\t' > reference/Mmus_ensembl.GO.csv
```

2. Map to species of interest using fine-grained orthologs:

```bash
Rscript s02_match_GOs_to_Mus_2023-08-22.R
```

#### Species tree of anthozoans

Obtain a species trees by selecting single-copy orthologs directly from Broccoli. This is simply used to assign branch lengths to the known species tree of Anthozoa. Thus, we are not terribly concerned about phylogenetic marker selection.

To run from the `data` folder.

1. All Anthozoa:

```bash
# species to include
spslist="Nvec Scocal Actieq Metsen Dialin Exapal Adig Amil Gfas Fspp Gasp Spis Ocupat Ocuarb Pocdam Xesp Dgig"
spsliststring="^Nvec\|^Scocal\|^Actieq\|^Metsen\|^Dialin\|^Exapal\|^Adig\|^Amil\|^Gfas\|^Fspp\|^Gasp\|^Spis\|^Ocupat\|^Ocuarb\|^Pocdam\|^Xesp\|^Dgig"
# dataset id
dcod="allanthozoa01"

# create output folder
mkdir -p orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments

# input data
for s in ${spslist} ; do
cat reference/${s}_long.pep.fasta 
done > orthology_Anthozoa_plus/phylogenomics_${dcod}/tmp.fasta
esl-sfetch --index orthology_Anthozoa_plus/phylogenomics_${dcod}/tmp.fasta

# get single-copy ogs
csvtool -t TAB -u TAB namedcol orthogroup,$(echo ${spslist} | tr ' ' ',') <(awk 'BEGIN {OFS="\t" } {if (NR==1) {$1="orthogroup"} } { print }' orthology_Anthozoa_plus/orthogroup_conservation.possvm.dollo_pres.csv) | awk '{for(i=2;i<=NF;i++){if($i == 1 && $(i-1)==$i){count++}};if((NF-2)==count){print};count=""}' | cut -f1 | head -n 200 > orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.txt
fgrep -f orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.txt <(awk 'BEGIN {OFS="\t" } {if (NR==1) {$1="orthogroup"} } { print }' orthology_Anthozoa_plus/orthogroup_conservation.possvm.dollo_pres.csv) > orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.csv

# get seqs and align
while read og ; do
ogs=$(echo $og | sed "s/://")
grep -w ${og} orthology_Anthozoa_plus/orthogroup_conservation.csv | cut -f4 | grep ${spsliststring} - > orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.txt
esl-sfetch -f orthology_Anthozoa_plus/phylogenomics_${dcod}/tmp.fasta orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.txt > orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.fasta
mafft --localpair --thread 8 --reorder --maxiterate 1000 orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.fasta > orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.l.fasta
clipkit orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.l.fasta -m kpic-gappy -o orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.lt.fasta -g 0.7
done < orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.txt

# concatenate
for s in ${spslist} ; do
echo ">${s}"
while read og ; do
ogs=$(echo $og | sed "s/://")
bioawk -c fastx '{ print $1,$2 }' orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.lt.fasta | grep "^${s}_" | head -n1 | cut -f2
done < orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.txt #| tr -d '\n'
echo ""
done > orthology_Anthozoa_plus/phylogenomics_${dcod}/dataset.fasta
bioawk -c fastx '{ print $1,length($2)} ' orthology_Anthozoa_plus/phylogenomics_${dcod}/dataset.fasta

# tree
iqtree -s orthology_Anthozoa_plus/phylogenomics_${dcod}/dataset.fasta -m TEST -mset LG -nt 4 -alrt 1000 -bb 1000 -pre orthology_Anthozoa_plus/phylogenomics_${dcod}/dataset.iqt -wsr -wbt
```

2. Species with available single-cell atlases:

```bash
# species to include
spslist="Ocupat Ocuarb Spis Amil Nvec Xesp"
spsliststring="^Ocupat\|^Ocuarb\|^Spis\|^Amil\|^Nvec\|^Xesp"
# dataset id
dcod="shortanthozoa01"

# create output folder
mkdir -p orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments

# input data
for s in ${spslist} ; do
cat reference/${s}_long.pep.fasta 
done > orthology_Anthozoa_plus/phylogenomics_${dcod}/tmp.fasta
esl-sfetch --index orthology_Anthozoa_plus/phylogenomics_${dcod}/tmp.fasta

# get single-copy ogs
csvtool -t TAB -u TAB namedcol orthogroup,$(echo ${spslist} | tr ' ' ',') <(awk 'BEGIN {OFS="\t" } {if (NR==1) {$1="orthogroup"} } { print }' orthology_Anthozoa_plus/orthogroup_conservation.possvm.dollo_pres.csv) | awk '{for(i=2;i<=NF;i++){if($i == 1 && $(i-1)==$i){count++}};if((NF-2)==count){print};count=""}' | cut -f1 | head -n 200 > orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.txt
fgrep -f orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.txt <(awk 'BEGIN {OFS="\t" } {if (NR==1) {$1="orthogroup"} } { print }' orthology_Anthozoa_plus/orthogroup_conservation.possvm.dollo_pres.csv) > orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.csv

# get seqs and align
while read og ; do
ogs=$(echo $og | sed "s/://")
grep -w ${og} orthology_Anthozoa_plus/orthogroup_conservation.csv | cut -f4 | grep ${spsliststring} - > orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.txt
esl-sfetch -f orthology_Anthozoa_plus/phylogenomics_${dcod}/tmp.fasta orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.txt > orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.fasta
mafft --localpair --thread 8 --reorder --maxiterate 1000 orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.fasta > orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.l.fasta
clipkit orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.l.fasta -m kpic-gappy -o orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.lt.fasta -g 0.7
done < orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.txt

# concatenate
for s in ${spslist} ; do
echo ">${s}"
while read og ; do
ogs=$(echo $og | sed "s/://")
bioawk -c fastx '{ print $1,$2 }' orthology_Anthozoa_plus/phylogenomics_${dcod}/alignments/${ogs}.lt.fasta | grep "^${s}_" | head -n1 | cut -f2
done < orthology_Anthozoa_plus/phylogenomics_${dcod}/single_copy_orthologs.txt #| tr -d '\n'
echo ""
done > orthology_Anthozoa_plus/phylogenomics_${dcod}/dataset.fasta
bioawk -c fastx '{ print $1,length($2)} ' orthology_Anthozoa_plus/phylogenomics_${dcod}/dataset.fasta

# tree
iqtree -s orthology_Anthozoa_plus/phylogenomics_${dcod}/dataset.fasta -m TEST -mset LG -nt 4 -alrt 1000 -bb 1000 -pre orthology_Anthozoa_plus/phylogenomics_${dcod}/dataset.iqt -wsr -wbt
```

#### Evolution of gene family gain/loss/expansion

From `results_annotation`:

1. Using Count:

```bash
mkdir -p results_gene_family_evolution_anthozoa

# get training dataset (from variable PFAM families)
Rscript s40_ancestral_reconstruction_prepare_dataset_2024-06-12.R

# Train with gamma categories for gain, loss, length and transfer (transfer seems important because it emulates capture errors)
# start without gamma cats
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -opt_delta 0.01 -max_paralogs 100 -uniform_duplication true -v true ../data/species_tree.Anthozoa.newick results_gene_family_evolution_anthozoa/pfam_domain_counts.train.csv > results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G0.txt
# add one gamma cat to the gain, loss, transfer, length and duplication parameters; and then iteratively expand to 2 and then 4
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -opt_delta 0.01 -max_paralogs 1000 -gain_k 1 -loss_k 1 -transfer_k 1 -length_k 1 -duplication_k 1 -uniform_duplication false -v true ../data/species_tree.Anthozoa.newick results_gene_family_evolution_anthozoa/pfam_domain_counts.train.csv results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G0.txt > results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G1.txt
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -opt_delta 0.01 -max_paralogs 1000 -gain_k 2 -loss_k 2 -transfer_k 2 -length_k 2 -duplication_k 2 -uniform_duplication false -v true ../data/species_tree.Anthozoa.newick results_gene_family_evolution_anthozoa/pfam_domain_counts.train.csv results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G1.txt > results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G2.txt
# java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -opt_delta 0.01 -max_paralogs 1000 -gain_k 4 -loss_k 4 -transfer_k 4 -length_k 4 -duplication_k 4 -uniform_duplication false -v true ../data/species_tree.Anthozoa.newick results_gene_family_evolution_anthozoa/pfam_domain_counts.train.csv results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G2.txt > results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G4.txt

## Same, but adding Ocuarb
# Train with gamma categories for gain, loss, length and transfer (transfer seems important because it emulates capture errors)
# start without gamma cats
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -opt_delta 0.01 -max_paralogs 100 -uniform_duplication true -v true ../data/species_tree.Anthozoa_plus.newick results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.train.csv > results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.rates.G0.txt
# add one gamma cat to the gain, loss, transfer, length and duplication parameters; and then iteratively expand to 2 and then 4
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -opt_delta 0.01 -max_paralogs 1000 -gain_k 1 -loss_k 1 -transfer_k 1 -length_k 1 -duplication_k 1 -uniform_duplication false -v true ../data/species_tree.Anthozoa_plus.newick results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.train.csv results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.rates.G0.txt > results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.rates.G1.txt
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -opt_delta 0.01 -max_paralogs 1000 -gain_k 2 -loss_k 2 -transfer_k 2 -length_k 2 -duplication_k 2 -uniform_duplication false -v true ../data/species_tree.Anthozoa_plus.newick results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.train.csv results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.rates.G1.txt > results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.rates.G2.txt
```

2. Run posterior probability analysis:

```bash
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.Posteriors -max_paralogs 10000 ../data/species_tree.Anthozoa.newick results_gene_family_evolution_anthozoa/orthogroup_counts.csv results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G0.txt  > results_gene_family_evolution_anthozoa/orthogroup_counts.posteriors.G0.csv
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.Posteriors -max_paralogs 10000 ../data/species_tree.Anthozoa.newick results_gene_family_evolution_anthozoa/orthogroup_counts.csv results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G1.txt  > results_gene_family_evolution_anthozoa/orthogroup_counts.posteriors.G1.csv
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.Posteriors -max_paralogs 10000 ../data/species_tree.Anthozoa.newick results_gene_family_evolution_anthozoa/orthogroup_counts.csv results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G2.txt  > results_gene_family_evolution_anthozoa/orthogroup_counts.posteriors.G2.csv
# java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.Posteriors -max_paralogs 10000 ../data/species_tree.Anthozoa.newick results_gene_family_evolution_anthozoa/orthogroup_counts.csv results_gene_family_evolution_anthozoa/pfam_domain_counts.rates.G4.txt  > results_gene_family_evolution_anthozoa/orthogroup_counts.posteriors.G4.csv

## Same, but adding Ocuarb
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.Posteriors -max_paralogs 10000 ../data/species_tree.Anthozoa_plus.newick results_gene_family_evolution_anthozoa_plus/orthogroup_counts.csv results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.rates.G0.txt  > results_gene_family_evolution_anthozoa_plus/orthogroup_counts.posteriors.G0.csv
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.Posteriors -max_paralogs 10000 ../data/species_tree.Anthozoa_plus.newick results_gene_family_evolution_anthozoa_plus/orthogroup_counts.csv results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.rates.G1.txt  > results_gene_family_evolution_anthozoa_plus/orthogroup_counts.posteriors.G1.csv
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.Posteriors -max_paralogs 10000 ../data/species_tree.Anthozoa_plus.newick results_gene_family_evolution_anthozoa_plus/orthogroup_counts.csv results_gene_family_evolution_anthozoa_plus/pfam_domain_counts.rates.G2.txt  > results_gene_family_evolution_anthozoa_plus/orthogroup_counts.posteriors.G2.csv
```

3. Parse gains/losses/expansions, and create tables of gene ages for each reference species:

```bash
Rscript s41_ancestral_reconstruction_posterior_analysis_2024-06-12.R
Rscript s43_ancestral_gain_loss_Dollo_2022-11-14.R # same, but with dollo (mostly consistent)
```

4. Parse gains and expansions per node, compute functional enrichments for specific nodes and along ancestral paths (note: excludes Hexacorallia node for scleractinians):

```bash
Rscript s42_ancestral_reconstruction_posterior_enrichments_2024-06-14.R
```

#### Analysis of tandem duplications

Analyse the presence of tandem duplications in extant anthozoan genomes.

To run from `results_annotation`.

1. Identify sets of extant tandem duplications based on ancestrally single-copy families at various ancestral nodes (Anthozoa, Hexacorallia, Scleractinia):

```bash
Rscript s44_segmental_duplications_from_ancestor_2024-07-04.R
```

2. Explore expression profile of segmentally duplicated families (only for species for which we have single-cell atlases; see below: [Single-cell transcriptome analysis](#single-cell-transcriptome-analysis)).

```bash
Rscript s45_segmental_duplications_expression_2024-07-05.R
```

#### Microsynteny analyses

Analyse patterns of microsynteny (gene pair collinearity) across anthozoan genomes.

To run from `results_annotation`.

1. Microsynteny conservation (pairs of syntenic genes):

```bash
Rscript s30_synteny_pairs_2024-05-27.R
Rscript s31_synteny_pairs_distance_2024-05-27.R
```

2. Synteny conservation in chromosomes (one-to-one ortholog placement across homolgous chromosomes):

```bash
Rscript s32_synteny_pairwise_2024-05-30.R
```

3. Explore expression profile of microsyntenic gene pairs (only for species for which we have single-cell atlases; see below: [Single-cell transcriptome analysis](#single-cell-transcriptome-analysis)).

```bash
## TODO
```

#### Macrosynteny analyses

Analyse patterns of macrosynteny, i.e. the evolution of ancestral linkage groups of genes (ALGs) in extant genomes.

To run from `results_annotation`.

1. Ancestral linkeage group analysis:

```bash
# 1. Align all-to-all, all species of interest:
# blast all species
mkdir -p results_macrosynteny_plus/alignments
mkdir -p results_macrosynteny_plus/bins
while read i ; do 
while read j ; do 
if [ $i != $j -a ! -f results_macrosynteny_plus/alignments/dmd.${i}-${j}.csv ] ; then 
echo $i $j
diamond blastp -d ../data/reference/${j}_long.pep.fasta -q ../data/reference/${i}_long.pep.fasta -o results_macrosynteny_plus/alignments/dmd.${i}-${j}.csv --evalue 1e-9 --more-sensitive --max-target-seqs 10 --quiet --threads 20
fi
done < data/species_list_synteny_blocks_plus.txt
done < data/species_list_synteny_blocks_plus.txt

# 2. Get clusters with MCL
while read i ; do 
while read j ; do 
if [ $i != $j ] ; then 
	awk '{print $1,$2 }' results_macrosynteny_plus/alignments/dmd.${i}-${j}.csv 
fi
done < data/species_list_synteny_blocks_plus.txt
done < data/species_list_synteny_blocks_plus.txt > results_macrosynteny_plus/mcl.all.abc.csv
mcl results_macrosynteny_plus/mcl.all.abc.csv --abc -I 2.1 -o results_macrosynteny_plus/mcl.all.out.csv
awk '{ { for(i = 1; i <= NF; i++) { printf("HG%06d\t%s\n", NR,$i) } } }' results_macrosynteny_plus/mcl.all.out.csv > results_macrosynteny_plus/mcl.all.out.txt

# 3. Match orthologroups to running windows of genes along each chromosome:
Rscript s60_define_regions_2022-05-30-whole.R

# 4. Define ancestral linkage groups based on homology to three reference species:
Rscript s61_find_ALG_2022-06-03.R

# 5. Score ancestral linkage groups in each species:
Rscript s62_score_ALG_2022-06-03.R

# 6. Match species chrs, pairwise
Rscript s63_pairwise_species_2022-06-08.R
```

#### Whole-genome alignment and conservation

Create whole-genome alignments using [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) and parse them to calculate Phastcons and PhyloP scores using [RPHAST](https://github.com/CshlSiepelLab/RPHAST).

To run from `results_annotation`.

1. Cactus alignments:

```bash
mkdir -p results_wga

# 1. Index, scleractinian corals with chromosomal genomes and oculina relatives
cat ../data/species_tree.Scleractinia_plus.newick > results_wga/index_cactus.Scleractinia_plus.txt
for i in Ocupat Ocuarb Fspp Gasp Pocver Spis Amil Acrpal Acrcer Porlut ; do
echo -e "${i}\t../data/reference/${i}_gDNA.fasta" >> results_wga/index_cactus.Scleractinia_plus.txt
done

# 2. Run cactus for scleractinian corals:
source ~/Programes/cactus-bin-v2.9.7/venv-cactus-v2.9.7/bin/activate
cactus results_wga/results_wga_cactus2 results_wga/index_cactus.Scleractinia_plus.txt results_wga/cactus_Scleractinia_plus.hal

# 3. Format conversion: hal to MAF
mkdir -p results_wga/results_alignments/
s="Ocuarb"
while read c ; do
hal2maf --refSequence ${c} --refGenome ${s} results_wga/cactus_Scleractinia_plus.hal results_wga/results_alignments_plus/cactus.${s}.${c}.maf
done < <(bioawk -c fastx 'length($2) > 1e6 { print $1 }' ../data/reference/${s}_gDNA.fasta)

```

2. Parse Cactus alignments:

```bash
# Pairs of homologous chromosomes
Rscript s70_cactus_pairs_of_scaffolds_2022-08-10.R

# Use MAF alignments to create neutral models for each species with `phyloFit`
Rscript s71_phastcons_2023-07-07.R

# Identify conservation and acceleration/conservation scores with `phastCons` and `phyloP`:
Rscript s72_phastcons_run_2023-07-18.R

# Identify conservation (presence) as percent identity across species, from MAF alignments
Rscript s73_cactus_parse_2023-07-18.R
```

#### Ks rate analysis of whole-genome duplications

Test for whole-genome duplications using paralog Ks distributions, using [`ksrates`](https://ksrates.readthedocs.io/). 

To run from `results_annotation`.

1. Launch `ksrates`:

```bash
mkdir -p results_wgd
conda activate ksrates
cd results_wgd
bioawk -c gff '$6=="."' ../../data/reference/Ocupat_long.annot.gtf | gffread - > Ocupat.gff3
bioawk -c gff '$6=="."' ../../data/reference/Amil_long.annot.gtf | gffread - > Amil.gff3
bioawk -c gff '$6=="."' ../../data/reference/Spis_long.annot.gtf | gffread - > Spis.gff3
for i in Nvec Scocal Actieq Metsen Dialin Exapal Adig Amil Gfas Fspp Gasp Spis Ocupat Pocdam Xesp Dgig ; do 
cp ../../data/reference/${i}_long.cds.fasta .
done
~/Programes/ksrates/nextflow run VIB-PSB/ksrates --config config_ksrates.Oculina.txt
~/Programes/ksrates/nextflow run VIB-PSB/ksrates --config config_ksrates.Acropora.txt
~/Programes/ksrates/nextflow run VIB-PSB/ksrates --config config_ksrates.Stylophora.txt
```

#### Transposon complement analysis

1. TE annotation using EDTA

```bash
## TODO
```

2. List types of elements and extract types per species, and align to self:

```bash
# alternative, aligning to consensus
mkdir -p results_repeat_analysis_plus/
for i in ../data/reference/*TE.all.gff3 ; do grep -v "^#" $i ; done | cut -f3 | sort | uniq -c | sort -nr | awk '{ print $2, $1 }' | tr ' ' '\t' > results_repeat_analysis_plus/repeat_counts.txt
cut -f1 results_repeat_analysis_plus/repeat_counts.txt | tr '_' '\t' | awk '{print $NF "\t" $0}' | sort -k1,1 | cut -f 2- | tr '\t' '_' > results_repeat_analysis_plus/repeat_list.txt

for i in Ocupat Ocuarb Fspp Gasp Spis Pocdam Pocver Adig Amil Acrpal Acrcer Gfas Porlut Actieq Metsen Dialin Exapal Nvec Scocal Dgig Xesp ; do 
rm -rf results_repeat_analysis_plus/te_complement_${i}
mkdir -p results_repeat_analysis_plus/te_complement_${i}
grep -v "^#" ../data/reference/${i}.TE.all.gff3 | cut -f3 | sort | uniq -c | sort -nr | awk '{ print $2, $1 }' | tr ' ' '\t' > results_repeat_analysis_plus/te_complement_${i}/repeat_counts.txt
makeblastdb -dbtype nucl -parse_seqids -in ../data/reference/EDTA_TE_libraries/renamed_genome.${i}.fasta.mod.EDTA.TElib.fa
while read ti ; do
echo "> process $i $ti..."
bioawk -c gff '$3 == "'${ti}'"' ../data/reference/${i}.TE.all.gff3 | bedtools getfasta -bed - -fi ../data/reference/${i}_gDNA.fasta | bioawk -c fastx '{print ">'${i}.${ti}.'"NR" "$1"\n"$2 }' > results_repeat_analysis_plus/te_complement_${i}/tes.${i}.${ti}.fasta
blastn -db ../data/reference/EDTA_TE_libraries/renamed_genome.${i}.fasta.mod.EDTA.TElib.fa -query results_repeat_analysis_plus/te_complement_${i}/tes.${i}.${ti}.fasta -out results_repeat_analysis_plus/te_complement_${i}/tes.${i}.${ti}.blast_to_fams.csv -num_threads 20  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qseq sseq"  -max_target_seqs 10
gzip results_repeat_analysis_plus/te_complement_${i}/tes.${i}.${ti}.blast_to_fams.csv
done < <(cut -f1 results_repeat_analysis_plus/te_complement_${i}/repeat_counts.txt)
done 
```

3. Calculate Kimura distance from alignments, plot distributions:

```bash
Rscript s51_kimura_from_tes_2024-11-20.R
```

### Single-cell transcriptomic atlas analyses

All the following commands and scripts are to be run from the `results_scatlas` folder, unless otherwise specified.

#### Mapping scRNAseq to corals

Using Cellranger `vXX`.

```bash
## TODO
```

#### Cluster and annotate scRNAseq atlases

The following scripts describe the process of cell clustering and cell type annotation for each species (`Ocupat`, `Spin`, `Amil`; and also the outgroups `Xesp`, `Ocuarb` and `Nvec`). For each species, results are organised a preliminary folder (denoted as `results_metacell_[SPSID]_prefilt/`) and a final, filtered folder (denoted as `results_metacell_[SPSID]_filt/`).

1. Initialise clustering solution with a first iteration (`*prefilt/`):

```bash
# normalise data for marker selection based on the SCT procedure
# also, discard doublets based on CT (and, for Ocupat, aided by doubletfinder)
Rscript s03_seurat_cells_normalise_2024-04-17-Ocupat.R    # for Ocupat, Amil and Spis
Rscript s03_seurat_cells_normalise_2024-04-17-Outgroups.R # for Xesp and Nvec
# create initial metacell solution and cluster
Rscript s04_seurat_cells_process_2024-04-17.R
```

2. Annotate and reorder clusters from the first iteration manually, and save an annotation file for this second iteration to `*filt/`.

3. Analyse the atlas after manual curation of clusters (markers per cell type, enrichments, cytotrace, etc); goes to `filt/` folders:

```bash
Rscript s06_cluster_postfilter_2024-04-18.R
```

4. Differential gene expression analyses between pairs of cell types (selected):

```bash
Rscript s07_dge_analyses_2024-01-29.R
```

5. TF combinations defining each cell type:

```bash
Rscript s13_TF_combinations_2024-05-09.R
```

6. Identify gene modules, annotate them and compare their contents across species, using `WGCNA` at the metacell level:

```bash
# create gene modules per species
Rscript s20_gene_modules_2024-04-29.R
# annotate them and run functional enrichment analyses
Rscript s21_gene_modules_annot_2024-04-29.R
# compare module contents across species, run trees
Rscript s22_gene_modules_comparisons_2024-06-05.R
# dedicated analyses for the modules in host cells (select manually form eigengenes plots)
Rscript s23_gene_modules_host_cells_2024-06-10.R
```

7. Assign transcriptomes from FACS-sorted symbiont-positive cells to cell types in the reference datasets.

```bash
Rscript s30_transfer_alga_positive_cells_symbionts_2024-06-26.R # based on MARS experiments
Rscript s31_metacell_symbiont_counts_2024-05-09.R # based on 10x data
```

8. Evolution of cell type-specific gene expression programmes:

```bash
# prepare possvm-style list of markers (OGs) per cell type
conda activate base
Rscript s50_clean_ct_evol_prepare_dollo_2024-10-23.R

# dollo for OG cell type markers
conda activate possvm
for i in results_ct_evolution/anc.*extant.tsv ; do
python ../scripts/possvm_reconstruction.py -tree results_ct_evolution/small_tree_plus.newick -ort <(cut -f1,2 ${i}) -out ${i%%.extant.tsv}.possvm
done
# dollo for all OGs (genome-level)
fgrep -w -f <(tr ',' '\n' < results_ct_evolution/small_tree_plus.newick | tr -d '(' | sed "s/).*//") ../data/orthology_Anthozoa_plus/orthogroup_conservation.csv | cut -f 1,3 | sed "1 i\orthogroup\tgene" > results_ct_evolution/small_tree_plus.orthogroups.csv
python ../scripts/possvm_reconstruction.py -tree results_ct_evolution/small_tree_plus.newick -ort results_ct_evolution/small_tree_plus.orthogroups.csv -out results_ct_evolution/small_tree_plus.orthogroups.possvm

# evolutionary reconstructions
Rscript s51_clean_ct_evol_process_dollo_2024-10-24.R
```

### Cross-species single-cell transcriptome analyses

Compare cell types across species in various ways.

#### Cell type trees across species

To run from `results_scatlas/` folder.

1. Cross-species similarity of cell types based on ICC best-gene pairs and weighted pearson correlation:

```bash
Rscript s08_csps_icc_calculation_2024-04-26.R # find best gene pairs between species using ICC
Rscript s09_csps_similarity_2024-04-29.R      # csps similarity plots and shared genes
```

2. Build cell type trees (i.e. hierarchical clustering of cell types) and identify node-specific markers.

```bash
Rscript s10_cell_type_trees_2024-06-21.R      # cell type trees (UPGMA, PCA-based...)
```

#### Use SAMap to align cell types across species

To run from `results_csps/` folder.

1. Create blast databases for SAMap:

```bash
# blast dbs
mkdir -p results_blast/dbs
for si in Ocupat Ocuarb Spis Amil Nvec Xesp ; do
bioawk -c fastx '{ print $1,$2 }' ../data/reference/${si}_long.pep.fasta | awk 'NR==FNR { l[$1]=$2;next} { for(i = 1; i <= NF; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' ../data/reference/${si}_all.annot.pep.gene2tx.txt - | awk '{print ">"$1"\n"$2 }'> results_blast/dbs/${si}.fasta
makeblastdb -dbtype prot -parse_seqids -in results_blast/dbs/${si}.fasta
done

# blast alignments (first is reference)
for si in Ocupat Ocuarb Spis Amil Nvec Xesp ; do
for sj in Ocupat Ocuarb Spis Amil Nvec Xesp ; do
bash qsub_blast.sh results_blast/dbs/${si}.fasta results_blast/dbs/${sj}.fasta  results_blast/samap.${si}-${sj}.blast.csv blastp
done
done

for si in Ocupat Ocuarb Spis Amil Nvec Xesp ; do
for sj in Ocupat Ocuarb Spis Amil Nvec Xesp ; do
mkdir -p results_samap/blasts/${si}${sj}
cp results_blast/samap.${si}-${sj}.blast.csv results_samap/blasts/${si}${sj}/${sj}_to_${si}.txt
cp results_blast/samap.${sj}-${si}.blast.csv results_samap/blasts/${si}${sj}/${si}_to_${sj}.txt
done
done

# modify blast outputs to fit Spio/Spin...
for si in Spin ; do
for sj in Ocupat Ocuarb Amil Nvec Xesp ; do
mkdir -p results_samap/blasts/${si}${sj}
cp results_blast/samap.Spis-${sj}.blast.csv results_samap/blasts/${si}${sj}/${sj}_to_${si}.txt
cp results_blast/samap.${sj}-Spis.blast.csv results_samap/blasts/${si}${sj}/${si}_to_${sj}.txt
done
done

for si in Ocupat Ocuarb Amil Nvec Xesp ; do
for sj in Spin ; do
mkdir -p results_samap/blasts/${si}${sj}
cp results_blast/samap.${si}-Spis.blast.csv results_samap/blasts/${si}${sj}/${sj}_to_${si}.txt
cp results_blast/samap.Spis-${si}.blast.csv results_samap/blasts/${si}${sj}/${si}_to_${sj}.txt
done
done
```

2. Prepare SAMap matrices from Seurat objects:

```bash
conda activate base
Rscript s30_objects_from_seurat_2024-07-02.R
```

3. Initialise single-species SAM objects:

```bash
conda activate samap
python s31_samap_prepare_SAM_2023-03-20.py
```

4. Cross-species analyses using SAMap:

```bash
# aligment, all-to-all (various combinations)
conda activate samap
python s35_samap_objects_all2all_2024-07-03.py
python s36_samap_postalignment_all2all_2024-07-03.py
conda activate base
Rscript s37_samap_plots_all2all_2024-07-03.R # heatmaps and pairwise igraphs

```

The end!

```python
				⣀⡀           
			⠘⢷⣤⣿⡇        
	⢰⡗  ⢠⡀⣠⡄ ⠈⣿   ⢀        
⠸⢶⣤⣄⢿⡇  ⠈⣿⠏   ⣿⡀  ⣴⠟⠁       
	⠙⠻⣿⣦⡀⢸⡏    ⢹⣇⣼⣏⣀⣀⣠⣤⡦    
⢰⣶⡄  ⠘⢿⣿⣾⣧    ⣼⣿⠟⠉⠉⠉⢉⡀     
⠘⠷⠶⢿⣿⡄   ⠙⠿⣿⣦⣄⡀⣼⣿⠃ ⠰⣦⣀⣸⡇     
⠈⣿⣷       ⢀⠈⠛⠿⣿⣿⡏   ⠈⣹⡟  ⣀⣤⣄ 
⢠⡶⠟⠻⣿⣧⡀⣰⠏    ⢸⣿⡇ ⣀⣠⣾⣯⣶⣶⣾⡏⠁  
	⠈⠻⣿⣿⣀     ⠸⣿⣷⣿⣿⣿⣯⣉⡉⠉⠙⢷⡄  
	⠈⠙⢿⣷⣦⣄ ⣠⣾⣿⠋   ⠈⣩⡿⠷⣤    
	⠶⠶⣶⡶⠟⠛⠛⢿⣿⡿⠁    ⣰⡟⠁      
    ⠹⠇   ⣼⣿⠃        ⠉⠁       
		 ⣿⣿               
~~~~~~~~^^^^^~~~~~~~~~~~~~~~~            
```
