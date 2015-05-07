DIR="/mnt/scratch0/dfilippo/aligned"
samfiles=("P_aeruginosa_PAO1" "SRR1294122" "SRR445718" "K562_cytosol_LID8465_TopHat_v2" "NA12878_S1")
# samfiles=("SRR1294122" "SRR445718" "P_aeruginosa_PAO1")
OUTPUT=all_files.sizes
rm $OUTPUT
for ((i=0; i < ${#samfiles[@]}; i++))
do
	echo $i
	FILE=${samfiles[$i]}
	SRC=$DIR/$FILE

	# ls -l $SRC.sam.[eolrhuif]* $SRC.sam.k\=* > $SRC.referee.sizes
	echo "# SAM" >> $OUTPUT
	du -b $SRC.sam >> $OUTPUT
	echo "# BAM" >> $OUTPUT
	du -b $SRC.bam >> $OUTPUT
	echo "# Referee" >> $OUTPUT
	du -b $SRC.sam.[eolrhuif]*lz $SRC.sam.k\=*lz >> $OUTPUT
	echo "# Quip" >> $OUTPUT
	du -b $SRC.sam.qp >> $OUTPUT
	echo "# Deez" >> $OUTPUT
	du -b $SRC.dz >> $OUTPUT
done

python python/parse_sizes.py $OUTPUT
