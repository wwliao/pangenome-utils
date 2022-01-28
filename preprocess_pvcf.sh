#!/bin/bash
# Usage: preprocess_pvcf.sh <pVCF file> <sample name>
VCF=$1
SAMPLE=$2
MAXSIZE=50
MEM="10G"
GRAPH=$(basename $VCF .vcf.gz)
REF="/gpfs/gibbs/pi/ycgh/HPP/HPRC_PLUS/GRCh38/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

bcftools view -a -s $SAMPLE -Ou $VCF \
    | bcftools norm -f ${REF} -c s -m -any -Ou \
    | bcftools view -e 'GT="ref" | GT~"\."' -f 'PASS,.' -Ou \
    | bcftools sort -m ${MEM} -T ${SAMPLE}.${GRAPH}_sort_tmp/ -Ou \
    | bcftools norm -d exact -Oz -o ${SAMPLE}.${GRAPH}.norm.vcf.gz \
    && bcftools index -t ${SAMPLE}.${GRAPH}.norm.vcf.gz \
    && bcftools view -e "STRLEN(REF)>${MAXSIZE} | STRLEN(ALT)>${MAXSIZE}" \
                     -r ${CHROMS} -Oz -o ${SAMPLE}.${GRAPH}.max${MAXSIZE}.chr1-22.vcf.gz \
                     ${SAMPLE}.${GRAPH}.norm.vcf.gz \
    && bcftools index -t ${SAMPLE}.${GRAPH}.max${MAXSIZE}.chr1-22.vcf.gz \
    && rm ${SAMPLE}.${GRAPH}.norm.vcf.gz*
