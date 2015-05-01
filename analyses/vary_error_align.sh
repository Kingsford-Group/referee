EXE=STAR
for FILE in SRR1294122 # SRR445718 # SRR1294122
do
	INPUT=/data/SeqComp/fastq/rnaseq/$FILE
	OUTPUT=/mnt/scratch0/dfilippo/aligned/$FILE\_all
	THREADS="--runThreadN 16"
	for Z in 1 0.9 0.8 0.7 0.6 0.5
	do
		echo "Error rate $Z"
		# GENOME="--genomeDir /mnt/scratch0/dfilippo/star_index"
		GENOME_HG19="--genomeDir /mnt/scratch0/dfilippo/star_idx/human/hg19_extra_dz/"
		OUTFILTER="--outFilterScoreMin 0 --outFilterMatchNmin 0 --outFilterScoreMinOverLread $Z --outFilterMatchNminOverLread $Z --outFilterMultimapNmax 20"
		OUTFORMAT="--outSAMunmapped Within --outSAMattributes MD --outFileNamePrefix $OUTPUT.z=$Z."
		echo "$EXE $GENOME_HG19 --readFilesIn $INPUT.fastq $THREADS $OUTFILTER $OUTFORMAT"
		time $EXE $GENOME_HG19 --readFilesIn $INPUT.fastq $THREADS $OUTFILTER $OUTFORMAT
	done
done
