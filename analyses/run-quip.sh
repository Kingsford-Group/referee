#HUMAN_GENOME_HG19=/mnt/scratch0/dfilippo/genomes/human/hg19/human_genome.fa
DIR=/mnt/scratch0/dfilippo/aligned
HUMAN_GENOME_HG19=/mnt/scratch0/dfilippo/genomes/deez/all_chromosomes_hg_19.fa
PAERUG_GENOME=/data/genomes/bacterial/Pseudomonas_aeruginosa_PAO1_uid57945/NC_002516.fna
#FILE=P_aueruginosa_PA01.strip.align
#FILE=SRR445718.strip.align
#FILE=SRR1294122.strip.align
#FILE=NA12878\_S1
#FILE=K562_cytosol_LID8465_TopHat_v2
#INDIR=/mnt/scratch0/dfilippo/aligned/stripped_sam/only_aligned
echo "Genome: $HUMAN_GENOME_HG19"
for FILE in SRR445718 SRR1294122 # NA12878\_S1 # K562_cytosol_LID8465_TopHat_v2 # NA12878\_S1 # SRR445718 SRR1294122
do
	echo "Compressing: $FILE"
	/usr/bin/time -v quip -v -r $HUMAN_GENOME_HG19 $DIR/$FILE.sam 2> $DIR/$FILE.quip.log
done

