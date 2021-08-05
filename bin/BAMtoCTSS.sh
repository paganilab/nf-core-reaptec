#!/usr/bin/env bash
Ref_chrom=$1 #'/local/home/ubuntu/Ref/homSap_v33/GRCh38.p13.genome.chrom.sizes'
Ref_pro1=$2 #'/local/home/ubuntu/Ref/FANTOM5_hg38/homSap_hg38_fair+new_CAGE_peaks_phase1and2.anno.bed'
Ref_enh=$3 #'/local/home/ubuntu/Ref/FANTOM5_hg38/homSap_hg38_permissive_enhancer.bed'

# mkdir ./Enhancer
# mkdir ./Promoter
BAM=$4
CPUS=$5
# for file1 in *.bam; do
NAME=$(basename $BAM .bam)

  # mkdir ./$NAME
  # mv $BAM ./$NAME
  # cd ./$NAME

C2=`echo "$NAME"_mq20_mappedcount.txt`

samtools view -@ $CPUS -bq 20 -cF 4 $BAM > $C2
 
FWBED=`echo "$NAME"_mq20.ctss.fwd.bedGraph`
REVBED=`echo "$NAME"_mq20.ctss.rev.bedGraph`
ALLBED=`echo "$NAME"_mq20.ctss.ALL.bedGraph`

samtools view -@ $CPUS -F 4 -u -q 20 $BAM | genomeCoverageBed -ibam /dev/stdin -5 -bg -strand + | sort -k1,1 -k2,2n > $FWBED
samtools view -@ $CPUS -F 4 -u -q 20 $BAM | genomeCoverageBed -ibam /dev/stdin -5 -bg -strand - | sort -k1,1 -k2,2n > $REVBED
samtools view -@ $CPUS -F 4 -u -q 20 $BAM | genomeCoverageBed  -ibam /dev/stdin -5 -bg | sort -k1,1 -k2,2n > $ALLBED

  samtools view -@ $CPUS $BAM -q 20 -F 0x104 -u \
    | bamToBed -i \
    | awk 'BEGIN{OFS="\t"}{
      if($6=="+") { print $1, $2,   $2+1, ".", $5, $6 }
      else        { print $1, $3-1, $3,   ".", $5, $6 }
      }' \
    | sort -k1,1 -k2,2n -k6,6 \
    | bedtools groupby -g 1,2,3,4,6 -c 1 -o count \
    | awk 'BEGIN{OFS="\t"}{tmp=$5; $5=$6; $6=tmp; print $0}' \
    | gzip -c > ${NAME}_mq20.ctss.bed.gz

 FWBW=`echo "$NAME"_mq20.fwd.bw`
 REVBW=`echo "$NAME"_mq20.rev.bw`
 ALLBW=`echo "$NAME"_mq20.all.bw`

 bedGraphToBigWig $FWBED $Ref_chrom $FWBW
 bedGraphToBigWig $REVBED $Ref_chrom $REVBW
 bedGraphToBigWig $ALLBED $Ref_chrom $ALLBW

 name=`echo "$NAME"_mq20`

  echo -e "01STAT:MAPPED"|paste /dev/stdin $C2 > .${name}.tmp

  cat $Ref_pro1 \
    | awk '{if($6 == "+"){print}}'\
    | bigWigAverageOverBed $FWBW /dev/stdin /dev/stdout \
    | cut -f 1,4 \
    >> .${name}.tmp

  cat $Ref_pro1 \
    | awk '{if($6 == "-"){print}}'\
    | bigWigAverageOverBed $REVBW /dev/stdin /dev/stdout \
    | cut -f 1,4 \
    >> .${name}.tmp

  sort .${name}.tmp > ${name}_promoter.count.txt

#cp ${name}.count.txt ../Promoter
 
#ENH=`echo "$NAME"_mq20_enhancer_count.txt`

cat $Ref_enh \
    | bigWigAverageOverBed $ALLBW /dev/stdin /dev/stdout \
    | cut -f 1,4 > ${name}_enhancer.count.txt

#cp $ENH ../Enhancer
  
#cd ./..

#done
