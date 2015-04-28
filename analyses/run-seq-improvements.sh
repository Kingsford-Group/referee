SETUP="0"

if [ $SETUP -eq 1 ]
then
	DIR=/mnt/scratch0/dfilippo/aligned
	DZDIR=/data/referee/deez_results/dz
	DEEZ=$DZDIR/deez
	AERUG_GENOME=/data/genomes/bacterial/Pseudomonas_aeruginosa_PAO1_uid57945/NC_002516.fna
	ECOLI_GENOME=/mnt/scratch0/dfilippo/genomes/bacterial/E_coli_DH10.fa
	# the same one use in Deez publication
	# HUMAN_GENOME1=/mnt/scratch0/dfilippo/genomes/human/hg19/human_genome.fa
	HUMAN_GENOME1=/mnt/scratch0/dfilippo/genomes/deez/all_chromosomes_hg_19.fa
	# "P_aeruginosa_PAO1" $AERUG_GENOME
	samfiles=("SRR1294122" "SRR445718" "K562_cytosol_LID8465_TopHat_v2" "MiSeq_Ecoli_DH10B_110721_PF")
	genomes=($HUMAN_GENOME1 $HUMAN_GENOME1 $HUMAN_GENOME1 $ECOLI_GENOME)
	TIME="-v"
else
	DIR=/data/referee/aligned
	DEEZ=/Users/giantlynx/deez/deez/deez
	AERUG_GENOME=$DIR/NC_002516_deez.fna
	HUMAN_GENOME1=$DIR/all_chromosomes_hg_19.fa
	samfiles=("NA12878_S1.sorted")
	genomes=($HUMAN_GENOME1)
	TIME="-lp"
fi



for ((i=0; i < ${#samfiles[@]}; i++))
do
    echo ""
    FILE=${samfiles[$i]}
    GENOME=${genomes[$i]}
    SRC=$DIR/$FILE.sam
    OUT=$DIR/$FILE
    LOG=$DIR/$FILE.seq-only.comp.log

    echo "----------------"
    echo " Processing $DIR/$FILE.sam, log at $DIR/$FILE.comp.log"
    echo " Genome: $GENOME"
    
    if [ ! -e $DIR/$FILE.sam ]
        then
        echo "File does not exist!"
        exit
    fi
    if [ ! -e $GENOME ]
        then
        echo "Genome does not exist!"
        exit
    fi

    STRIPPED=$DIR/$FILE.stripped
    # only strip if the file does not already exist
    if [ ! -e $STRIPPED.sam ]
        then
        echo " stripping file..."
        time ./strip-fields-from-sam.sh $SRC > $STRIPPED.sam
    fi
    # wc -l $DIR/$FILE.stripped.sam

    # run deez
    DZ_OUTPUT=$STRIPPED.dz
    rm -f $DZ_OUTPUT
    echo " Compressing w/ Deez: $FILE"
    /usr/bin/time $TIME $DEEZ -r $GENOME -t 8 $STRIPPED.sam -o $DZ_OUTPUT 2> $STRIPPED.deez.log
    ls -l $DZ_OUTPUT
    # decompress
    # /usr/bin/time $TIME $DEEZ -! -r $GENOME $OUTPUT -o $STRIPPED.dz_recovered.sam
    # ls -l $STRIPPED.dz_recovered.sam

    # run CRAM
    # CRAM_OUTPUT=$STRIPPED.cram
    # rm -f $CRAM_OUTPUT
    # echo " Compressing w/ CRAM: $FILE"
    # TODO

    # run quip
    QP_OUTPUT=$STRIPPED.sam.qp
    rm -f $QP_OUTPUT
    echo " Compressing w/ Quip: $FILE"
    /usr/bin/time $TIME quip -v -r $GENOME $STRIPPED.sam 2> $STRIPPED.quip.log
    ls -l $STRIPPED.sam.qp
done
