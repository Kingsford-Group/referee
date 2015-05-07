DZDIR=/data/referee/deez_results/dz
HUMAN_GENOME=/mnt/scratch0/dfilippo/genomes/deez/all_chromosomes_hg_19.fa
INDIR=/mnt/scratch0/dfilippo/aligned

echo "Genome: $HUMAN_GENOME"

for FILE in NA12878\_S1 # SRR445718 SRR1294122 K562_cytosol_LID8465_TopHat_v2 # NA12878\_S1 K562_cytosol_LID8465_TopHat_v2
do
	INPUT=$INDIR/$FILE.sam
	OUTPUT=$INDIR/$FILE.dz
	rm -f $OUTPUT
	echo "Compressing: $FILE"
	echo "/usr/bin/time -v $DZDIR/deez -r $HUMAN_GENOME $INPUT -o $OUTPUT 2> $INDIR/$FILE.dz_comp.log"
	/usr/bin/time -v $DZDIR/deez -t 10 -r $HUMAN_GENOME $INPUT -o $OUTPUT 2> $INDIR/$FILE.dz_comp.log
done
