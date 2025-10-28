#!/bin/bash
#
#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH -J kite_ADT
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=8 # Request that ncpus be allocated per process.
#SBATCH --mail-user=zhongli@pitt.edu #send email to this address if ...

fastq=/ihome/wchen/zhongli/link2rduerr/shared/rduerr_wchen/Single_cell_fastq/Duerr_20220429_DOGMAseq-2_fastq/*ADT*
cellranger=/ihome/wchen/zhongli/RWorkSpace/DOGMA-seq/Duerr_20220429_DOGMAseq-2_arc
cut -f 1 -d , $cellranger/per_barcode_metrics.csv | sed 's/-1//g' | grep barcode -v > barcodes.txt
kallisto bus -i ../FeaturesMismatch.idx -o ./ -x 10xv3 -t 8 $fastq
bustools correct -w ./barcodes.txt ./output.bus -o ./output_corrected.bus
bustools sort -t 32 -o ./output_sorted.bus ./output_corrected.bus
mkdir ./featurecounts/
bustools count -o ./featurecounts/featurecounts --genecounts -g ../FeaturesMismatch.t2g -e ./matrix.ec -t ./transcripts.txt ./output_sorted.bus
