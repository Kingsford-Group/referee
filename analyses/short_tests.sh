make clean
make
BIN=bin/referee
if [ ! -e $BIN ]; then
        echo "Binary not found: $EXE"
    exit
fi
DIR=test
#REF=/data/referee/aligned/human\_chr10.fa
REF=/data/genomes/Homo_sapiens/NCBI/build37.2/Sequence/WholeGenomeFasta/genome.fa
for FILE in ordered addtnl\_fields insertion\_test no\_md\_string
do
	echo '--------'
	echo $DIR/$FILE.sam
	./$BIN -c -i $DIR/$FILE.sam -o $DIR/out/$FILE.fast --fast 2> $DIR/$FILE.comp.log
	./$BIN -d -i $DIR/out/$FILE.fast -r $REF -o $DIR/$FILE.recovered.sam 2> $DIR/$FILE.decomp.log
	# echo "diff --suppress-common-lines -w -y -a -B <(cut -f 1-5,7-10 $DIR/$FILE.recovered.sam) <(cut -f 1-5,7-10 $DIR/$FILE.sam)"
	diff --suppress-common-lines -w -y -a -B <(cut -f 1-5,7-10 $DIR/$FILE.recovered.sam) <(cut -f 1-5,7-10 $DIR/$FILE.sam) | wc -l
done

# consistent\_edits


# test decompression of insertions
echo "----------------"
echo " HARD CLIPS TEST "
FILE=hard_clips
REF=/data/referee/aligned/NC_002516_deez.fna
$BIN -c -i $DIR/$FILE.sam -o $DIR/out/$FILE --fast 2> $DIR/$FILE.comp.log

$BIN -d -i $DIR/out/$FILE -r $REF -o $DIR/$FILE.recovered.sam 2> $DIR/$FILE.decomp.log
diff --suppress-common-lines -w -y -a -B <(cut -f 1-5,7-10 $DIR/$FILE.recovered.sam) <(cut -f 1-5,7-10 $DIR/$FILE.sam) | wc -l
