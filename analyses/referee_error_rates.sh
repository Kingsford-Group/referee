cd ..
BIN=bin/referee
./so_plzip.sh
make clean
make

if [ ! -e $BIN ]
	then
	echo "Binary not found"
	exit
fi

export LD_LIBRARY_PATH=plzip:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=plzip:$DYLD_LIBRARY_PATH

# sudo purge -n
TIMEOPT="-v"

OUTPUT=error_rates.txt
rm $OUTPUT # write anew each time
DATADIR=/mnt/scratch0/dfilippo/aligned
for FILE in SRR445718 SRR1294122
do
	for Z in 0.9 0.8 0.7 0.6 0.5
	do
		#INPUT=$DATADIR/$FILE\_all.z=$Z
		#ls -l $INPUT.Aligned.out.sam.sorted.sam
		#mv $INPUT.Aligned.out.sam.sorted.sam $INPUT.sam
		INPUT=$DATADIR/$FILE\_all.z=$Z.sam
		echo "-----------------"
		echo "Processing $INPUT"
		ls -l $INPUT
		/usr/bin/time $TIMEOPT $BIN -c $INPUT 10 --seq2 > $INPUT.referee.seq2.log 2>&1
		echo "*** $FILE z=$Z" >> $OUTPUT
	        du -b $INPUT.offs.lz >> $OUTPUT
	        du -b $INPUT.left_clip.lz >> $OUTPUT
	        du -b $INPUT.right_clip.lz >> $OUTPUT
	        du -b $INPUT.edits.lz >> $OUTPUT
		# du -b $INPUT.has_edits >> $OUTPUT
		plzip -n 10 $INPUT.has_edits
	        du -b $INPUT.has_edits.lz >> $OUTPUT
	done
done

echo "Aggregated unique_seq_only sizes"
python python/parse_compr_rates_error_rate.py $OUTPUT

cd analyses
