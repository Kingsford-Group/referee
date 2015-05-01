DIR=/mnt/scratch0/dfilippo/aligned
HUMAN_GENOME_HG19=/mnt/scratch0/dfilippo/genomes/deez/all_chromosomes_hg_19.fa
PAERUG_GENOME=/data/genomes/bacterial/Pseudomonas_aeruginosa_PAO1_uid57945/NC_002516.fna
#FILE=NA12878\_S1
#FILE=K562_cytosol_LID8465_TopHat_v2

samfiles=("P_aeruginosa_PAO1" "K562_cytosol_LID8465_TopHat_v2")
genomes=($PAERUG_GENOME $HUMAN_GENOME_HG19 $HUMAN_GENOME_HG19)
for ((i=0; i < ${#samfiles[@]}; i++))
do
	echo ""
	FILE=${samfiles[$i]}
	GENOME=${genomes[$i]}
	echo "Compressing: $FILE"
	rm $DIR/$FILE.sam.qp
	/usr/bin/time -v quip -v -r $HUMAN_GENOME_HG19 $DIR/$FILE.sam 2> $DIR/$FILE.quip.log
done

