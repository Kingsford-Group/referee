EXE=bin/referee
make clean
make
if [ ! -e $EXE ]; then
        echo "Binary not found: $EXE"
    exit
fi

./so_plzip.sh

export LD_LIBRARY_PATH=plzip:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=plzip:$DYLD_LIBRARY_PATH

# purge memory cache
# sudo purge -n

#TIMEOPT="-lp"
TIMEOPT="-v"

DIR="/mnt/scratch0/dfilippo/aligned"
# DIR="/data/referee/aligned"
# SEQMODE="--seqOnly" # --discardSecondary

# samfiles=("P_aeruginosa_PAO1" "K562_cytosol_LID8465_TopHat_v2" "NA12878_S1")
samfiles=("NA12878_S1")
# samfiles=("SRR1294122" "SRR445718" "P_aeruginosa_PAO1")
# samfiles=("SRR445718" "SRR1294122" "P_aeruginosa_PAO1" "K562_cytosol_LID8465_TopHat_v2" "NA12878_S1")
#samfiles=("SRR1294122")
#samfiles=("MiSeq_Ecoli_DH10B_110721_PF")
#samfiles=("P_aeruginosa_PAO1.10mil")
#GENOME=NC_002516.fna
HUMAN=/mnt/scratch0/dfilippo/genomes/deez/all_chromosomes_hg_19.fa
PAERUG=/data/genomes/bacterial/Pseudomonas_aeruginosa_PAO1_uid57945/NC_002516.fna
#genomes=($PAERUG $HUMAN $HUMAN $HUMAN $HUMAN)
genomes=($HUMAN)
for ((i=0; i < ${#samfiles[@]}; i++))
do
	echo $i
	#if [ $i -gt 0 ] # change to 1 if running w/ P aeruginosa
	#then
	#	DIR=${dir[1]}
	#else
	#	DIR=${dir[0]}
		# continue
	#fi

	FILE=${samfiles[$i]}
	GENOME=${genomes[$i]}
	echo $FILE, $GENOME

	# ls -l $DIR/$FILE.sam.[oefilrhmk]*.lz
	rm -f $DIR/$FILE.sam.[oefilrhmk]*.lz

	SRC=$DIR/$FILE.sam
	rm -f $DIR/$FILE.sam.k\=[0-9].*
	OUT=$DIR/$FILE
	#LOG=$DIR/$FILE.comp.log
	if [ -z "$SEQMODE" ]
		then
		LOG=$DIR/$FILE.comp.log
	else
		LOG=$DIR/$FILE.seq_comp.log
	fi

	# ls -l $SRC
	echo "----------------"
	echo " Compress $DIR/$FILE.sam, log at $LOG"
	/usr/bin/time $TIMEOPT $EXE -c $SRC -t 10 -r $GENOME $SEQMODE 2> $LOG
	plzip -fv -n 10 $SRC.has_edits
	plzip -fv -n 10 $SRC.*membership
done
