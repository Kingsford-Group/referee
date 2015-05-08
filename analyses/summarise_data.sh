DIR=/mnt/scratch0/dfilippo/aligned
BIN=./bin/rsupport
OUTPUT=data_summary.log
rm -f $OUTPUT
for file in SRR445718 SRR1294122 P_aeruginosa_PAO1 K562_cytosol_LID8465_TopHat_v2 NA12878_S1
do
	echo $file
	plzip -df -n 10 $DIR/$file.sam.edits.lz
	plzip -df -n 10 $DIR/$file.sam.has_edits.lz
	echo "=== $file" >> $OUTPUT
	$BIN edits $DIR/$file.sam 100 2>> $OUTPUT
	# compress these files again to leave everything as we found it
	plzip -n 10 $DIR/$file.sam.edits
	plzip -n 10 $DIR/$file.sam.has_edits

	# count unaligned
	UNALIGN=$DIR/$file.sam.unaligned
	plzip -fd -n 10 $UNALIGN.lz
	wc -l $UNALIGN >> $OUTPUT
	# compress again
	plzip -n 10 $UNALIGN
done

# python python/parse_dataset_summary.py
