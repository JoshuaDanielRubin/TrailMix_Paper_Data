#!/bin/bash
# Ensure the script exits on any error
set -e

# Define the input file and base name
input_file=$1
base=$(basename $input_file .bam)

# Index the input BAM file
samtools index ${input_file}

# Get the mtDNA mapped reads (dummy step, use actual extraction if needed)
cp ${input_file} ${base}_mt.bam
echo "########### mtDNA bam OK ###########"

# Generate VCF file with variable sites
samtools mpileup -d 2000 -m 3 -C50 -q 30 -EQ 20 -uf ../human_MT.fa ${base}_mt.bam | bcftools call -m --ploidy 1 > ${base}_mt.vcf

echo "########### vcf OK ###########"

# Generate consensus sequence
./CnsMaj3_1.pl -i ${base}_mt.vcf -o ${base}.fa -l 16569 -cov 1 -diff 0.001 -idiff 0.001 -h ${base} -callindels no > ${base}_not_used.fa
echo "########### consensus file OK ###########"


# Align the consensus with 311 mtDNA sequences
cat ${base}.fa 311humans.fasta > ${base}_311humans.fasta
/home/ctools/mafft-7.487-with-extensions/scripts/mafft ${base}_311humans.fasta > aligned.${base}_311humans.fasta
echo "########### alignment with 311 mtDNAs OK ###########"

# Index the consensus sequence
bwa index -a bwtsw ${base}.fa
samtools faidx ${base}.fa
java -jar /home/ctools/picard_2.27.5/picard.jar CreateSequenceDictionary R=${base}.fa O=${base}.dict
echo "########### consensus indexed OK ###########"

# Generate FASTQ from BAM
java -jar /home/ctools/picard_2.27.5/picard.jar SamToFastq INPUT=${base}_mt.bam FASTQ=${base}.fq
#samtools bam2fq ${base}_mt.bam > ${base}.fq

echo "########### fastq from bam OK ###########"

# Map the FASTQ file to the consensus
bwa mem -R "@RG\tID:${base}\tLB:${base}_L1\tPL:ILLUMINA\tSM:${base}" ${base}.fa ${base}.fq | samtools view -bh -q 30 | samtools sort -O BAM -o ${base}_ra.sort.bam
echo "########### map fastq to the consensus OK ###########"

# Remove duplicates and re-index
java -jar /home/ctools/picard_2.27.5/picard.jar MarkDuplicates I=${base}_ra.sort.bam O=${base}_ra.sort.rmdup.bam TMP_DIR=temp METRICS_FILE=${base}.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
samtools index ${base}_ra.sort.rmdup.bam

# Use the consensus as the reference for calmd
samtools calmd -Erb ${base}_ra.sort.rmdup.bam ${base}.fa > ${base}_ra.final.bam
samtools index ${base}_ra.final.bam
echo "########### calmd step completed ###########"

# Run the contamMix analysis
../exec/estimate.R --samFn ${base}_ra.final.bam --malnFn aligned.${base}_311humans.fasta --figure ${base}_contam.pdf --trimBases 2 --nIter 100000 --nChains 5  | tee ${base}.summary.txt
