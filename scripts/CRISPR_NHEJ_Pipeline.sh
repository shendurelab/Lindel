#!/bin/bash
# Load parameters
samplename=$1
UMIs=$2
ref=$3
outdir=$4${samplename}
read_fwd=$5
read_rev=$6
bc_fwd=$7
bc_rev=$8
umi_correction=$9 
threshold=${10}
merge=${11}
workdir=${12}

# UMI_extract and process
echo "Runing UMI_extractor" 
python ${workdir}scripts/UMI_extractor.py -S $samplename -F $read_fwd -R $read_rev -O $outdir --bc_fwd $bc_fwd --bc_rev $bc_rev --umi_correction $umi_correction --threshold $threshold

# Merge and align two fastq files
if [ "$merge" == "TRUE" ]; then
	echo "Merging paired end reads."
	#pear-0.9.6-bin-32 -f $outdir/${samplename}_fwd.fastq -r $outdir/${samplename}_rev.fastq -o $outdir/${samplename}
	#rm $outdir/${samplename}.discarded.fastq $outdir/${samplename}.unassembled.forward.fastq $outdir/${samplename}.unassembled.reverse.fastq 
	#needleall $ref $outdir/${samplename}.assembled.fastq -datafile /net/shendure/vol10/projects/CRISPR.OT.Homing/nobackup/ref/EDNAFULL.modified -gapopen 20 -gapextend 0.5 -aformat3 fasta -awidth3=5000 -outfile $outdir/${samplename}_needleall.fasta
	python ${workdir}scripts/NHEJ_analysis.py $outdir ${samplename} 254
fi