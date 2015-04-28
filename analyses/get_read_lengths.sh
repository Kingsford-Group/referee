# generate read length and # alignments table
# DIR=/data/referee/aligned
DIR=/mnt/scratch0/dfilippo/aligned
OUT=read_lengths.all
rm -f $OUT
for FILE in SRR445718.stripped SRR1294122.stripped MiSeq_Ecoli_DH10B_110721_PF.stripped P_aeruginosa_PAO1.stripped K562_cytosol_LID8465_TopHat_v2.stripped NA12878_S1.stripped
do
	if [ ! -e $DIR/$FILE.sam ]
		then
		echo "File does not exist: $DIR/$FILE.sam" 
		continue
	fi
	echo $FILE >> $OUT;
	head -100 $DIR/$FILE.sam | tail -1 | cut -f 10 | wc -c >> $OUT
	wc -l $DIR/$FILE.sam | awk '{print $1}' >> $OUT
done
echo "Output written to $OUT"
