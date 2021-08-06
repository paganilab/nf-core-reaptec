#!/usr/bin/env bash

# Developed by Akiko Oguchi <aoguchi@kuhp.kyoto-u.ac.jp>
# Adapted by R.Bonnal <raoul.bonnal@ifom.eu> to work with Nextflow
# removing the for loop on bam files

## How to run this Script, for example :
## bash thisScript.sh G > xxx.log 2>&1
## or
## bash thisScript.sh A > xxx.log 2>&1


BAM=$1
BASE=$2
CPUS=$3

if   test ${BASE} == "A"; then
    ALT=T
elif test ${BASE} == "G"; then
    ALT=C
elif test ${BASE} == "C"; then
    ALT=G
elif test ${BASE} == "T"; then
    ALT=A
fi

NAME=`echo "${BAM}"|cut -d'.' -f1`

samtools view -H ${BAM} > ${NAME}_header.sam
#Gain reads with one-nucleotide soft-clipped on the forward strand
samtools view -@ ${CPUS} -F 16 ${BAM} | awk -F '\t' -v BASE=${BASE} '
    BEGIN {OFS="\t"} {
      if ($6 ~ /^1S[0-9]/ && $10 ~ /^'${BASE}'/) {print $0}
    }' \
    > ${NAME}_Softclip${BASE}_F.sam
#Gain reads with one-nucleotide soft-clipped on the reverse strand
samtools view -@ ${CPUS} -f 16 ${BAM} | awk -F '\t' -v ALT=${ALT} '
    BEGIN {OFS="\t"} {
      if ($6 ~ /[0-9]M1S$/ && $10 ~ /'${ALT}'$/) {print $0}
    }' \
    > ${NAME}_Softclip${BASE}_R.sam

#Unite the header, F.sam and R.sam
cat \
    ${NAME}_header.sam \
    ${NAME}_Softclip${BASE}_F.sam \
    ${NAME}_Softclip${BASE}_R.sam \
    | samtools sort -@ ${CPUS} -O bam -o Softclip${BASE}_${NAME}.bam

#Delete unnecessary files
rm ${NAME}*.sam
