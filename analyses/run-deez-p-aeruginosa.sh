DIR=/data/referee/aligned
GENOME=$DIR/NC_002516_deez.fna
FILE=P_aeruginosa_PAO1.100k.stripped
SRC=$DIR/$FILE.sam
OUTPUT=$DIR/$FILE.dz
DEEZ=/Users/giantlynx/deez/deez/deez
rm -f $OUTPUT
echo "Genome: $GENOME"
echo "Compressing: $FILE"
ls -l $SRC
/usr/bin/time -lp $DEEZ -r $GENOME -t 8 $SRC -o $OUTPUT
ls -l $OUTPUT
# decompress
/usr/bin/time -lp $DEEZ -r $GENOME $OUTPUT -o $DIR/$FILE.dz_recovered.sam
ls -l $DIR/$FILE.dz_recovered.sam