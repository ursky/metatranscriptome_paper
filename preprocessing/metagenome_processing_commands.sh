#!/usr/bin/env bash
# These are the commands used for processing the whole metagenomic sequencing data from the halite multi-omic study.
# The pipeline assumes that metaWRAP v1.1 is installed with default settings, and is available in the PATH.
# It is assumed that the 12 fastq raw read files for the metagenomic (DNA) sequencing of the samples are are in the "RAW_READS" directory.

# 1. Preprocessing the reads, trimming, human read contaminaiton removal:
tar -xvf RAW_READS/*tar.gz
mkdir READ_QC
for i in RAW_READS/*_1.fastq; do
	metawrap read_qc -1 ${i%_*}_1.fastq -2 ${i%_*}_2.fastq -t 12 -o READ_QC/${i%_*}
done

mkdir CLEAN_READS
for i in READ_QC/*; do
	b=${i#*/}
	mv ${i}/final_pure_reads_1.fastq CLEAN_READS/${b}_1.fastq
	mv ${i}/final_pure_reads_2.fastq CLEAN_READS/${b}_2.fastq
done

# At this point make sure that the files in "CLEAN_READS" are named as something meaningfull. For this example, I will rename the samples with the naming convention 9AM-1, 9AM-2, etc.


# 2. Metagenomic co-assembly of all the samples with metaSPAdes:
cat CLEAN_READS/*_1.fastq >> CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/*_2.fastq >> CLEAN_READS/ALL_READS_2.fastq
metawrap assembly --metaspades -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -m 1000 -t 112 -o ASSEMBLY
# With 112 threads, this took ~24 hours. Signifficant memory allocaiton (500BG+) was also necessary. If this is not possible, use MegaHIT instead (--megahit option), which requires signifficantly less resources and produces a comparable assembly.


# 3. Binning of scaffolds to produce metagenome assembled genomes:
# initial binning with metabat2, concoct, and maxbin2:
metawrap binning -o INITIAL_BINNING -t 112 -a ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct CLEAN_READS/*fastq
# inerative refinement with metaWRAP bin_refinement module (-c defines minimum desided completion, and -x maximum contamination):
metawrap bin_refinement -o BIN_REFINEMENT -t 112 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 70 -x 5
# MAG reassembly with the Reassemble_bins module was skipped for this analysis


# 4. MAG and contig abundance estimation. 
# This produces both the MAG abundance matrix (QUANT_BINS/abundance_table.tab) and the contig abundance matrix (QUANT_BINS/contig_abundance.tab).
metawrap quant_bins -b BIN_REFINEMENT/metaWRAP_bins -o QUANT_BINS -a ASSEMBLY/final_assembly.fa CLEAN_READS/9*fastq


# 5. Functional annotation:
# The annotation of the entire metagenomic assembly and the MAGs were obtained through the Integrated Microbial Genomes and Microbiomes (IMG) annotation service (through JGI). See manuscript for data availability and hosting details.
# In-house scripts (found under metatranscriptome_paper/funcitonal_annotation_processing) was used to rename the contigs in the IMG product_names file to the original contig names, and then produce a final img_annotation.master file containing the annotation and metadata for each gene in the assembly.


# 6. Taxonomic classification of contigs:
# The taxonomic predictions for the annoted genes produced from the JGI IMG annotation service (see img_annotation.taxonomy) was used to classify the contigs. A in-house script (combine_gene_taxonomy.py) was used to consolidate the gene taxonomies to predict the most likely annotation for the entire contig.


# 7. MAG classification:
# First pass MAG taxonomic classification
metawrap classify_bins -b BIN_REFINEMENT/metaWRAP_bins -o CLASSIFICATION -t 48
# More accurate classifications of the MAGs was obtained by manually aligning annotated sequences from each bin to the RefSeq database. 
# Sequences used were: 16S, 18S, 23S, 28S (nucleotide sequences), and ribosomal proteins L1, L2, L3, L4, S1, S2, S3, S4 (amino acid sequences). 
# Additional information was gained by manually loogin over the taxonomy assignemnts of the contigs belonging to each MAG. 

# 8. Metatranscriptome read processing:
# Trim and curate metatranscriptome reads exactly like the metagenomic reads in step #1

# 9. Quantify MAG and contig total transcription, and gene expression:
# This produces both the MAG total activity matrix (QUANT_BINS/abundance_table.tab) and the contig activity matrix (QUANT_BINS/contig_abundance.tab).
# For the gene expression estimation, use the transcript fasta file produced from the metagenome assembly functional annotation.
metawrap quant_bins -b BIN_REFINEMENT/metaWRAP_bins -o MAG_ACTIVITY -a ASSEMBLY/final_assembly.fa RNA_READS/9*fastq
metawrap quant_bins -b BIN_REFINEMENT/metaWRAP_bins -o GENE_ACTIVITY -a img_annotation.fa RNA_READS/9*fastq




