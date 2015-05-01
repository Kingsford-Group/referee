EXE=bin/referee
if [ ! -e $EXE ]; then
        echo "Binary not found: $EXE"
    exit
fi

./so_plzip.sh

export LD_LIBRARY_PATH=plzip:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=plzip:$DYLD_LIBRARY_PATH

# purge memory cache
# sudo purge -n

SEQMODE="--seq"

DIR="/mnt/scratch0/dfilippo/aligned"
# samfiles=("P_aeruginosa_PAO1" "SRR1294122" "SRR445718" "K562_cytosol_LID8465_TopHat_v2")
# samfiles=("SRR1294122" "SRR445718" "P_aeruginosa_PAO1")
#samfiles=("SRR1294122" "SRR445718" "K562_cytosol_LID8465_TopHat_v2")
#samfiles=("SRR1294122")
#samfiles=("MiSeq_Ecoli_DH10B_110721_PF")
samfiles=("P_aeruginosa_PAO1")
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
	/usr/bin/time -v $EXE -c $SRC 10 $SEQMODE 2> $LOG
	plzip -f -n 10 $SRC.has_edits
	#plzip -f -n 10 $SRC.*membership
done
