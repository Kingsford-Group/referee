BIN=./bin/referee
# ./so_plzip.sh
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

# DATADIR=/mnt/scratch0/dfilippo/aligned
DATADIR=/data/referee/aligned
# FILE=SRR445718.sam
FILE=SRR445718.sam
#FILE=SRR1294122.sam
#DATADIR=.
# FILE=test.sam
#FILE=P_aeruginosa_PAO1.10mil.sam
#FILE=K562_cytosol_LID8465_TopHat_v2.sam
# FILE=NA12878_S1.sorted.sam
rm -f $DATADIR/$FILE.offs.lz
rm -f $DATADIR/$FILE.edits.lz
rm -f $DATADIR/$FILE.left_clip.lz
rm -f $DATADIR/$FILE.right_clip.lz
rm -f $DATADIR/$FILE.flags.lz
rm -f $DATADIR/$FILE.ids.lz
rm -f $DATADIR/$FILE.quals.lz
rm -f $DATADIR/$FILE.opt.lz
rm -f $DATADIR/$FILE.unaligned.lz
rm -f $DATADIR/$FILE.k\=*

# TIMEOPT="-v"
TIMEOPT="-lp"

/usr/bin/time $TIMEOPT $BIN -c $DATADIR/$FILE 10 --seq2 > $DATADIR/$FILE.referee.log 2>&1
# ls -l $DATADIR/$FILE.*
# time plzip -vf $DATADIR/$FILE.k\=* 2> $DATADIR/$FILE.plzip
# python python/parse_plzip_output.py $DATADIR/$FILE.plzip
# $BIN -c $DATADIR/$FILE 10
# ls -l $DATADIR/$FILE.*
# /usr/bin/time -lp $BIN -c $DATADIR/P_aeruginosa_PAO1.10mil.sam

# DECOMPRESS=1
# if [ $DECOMPRESS -eq 1 ]
# then
# 	HUMAN_GENOME=$DATADIR/human_genome.fa
# 	# HUMAN_GENOME=$DATADIR/../genomes/deez/all_chromosomes_hg_19.fa
# 	plzip -fd $DATADIR/$FILE.offs.lz
# 	plzip -fd $DATADIR/$FILE.edits.lz
# 	plzip -fd $DATADIR/$FILE.*clip.lz
# 	plzip -fd $DATADIR/$FILE.has_edits.lz
# 	/usr/bin/time $TIMEOPT $BIN -d $DATADIR/$FILE $HUMAN_GENOME > $DATADIR/$FILE.referee.decomp.log 2>&1
# 	# diff --suppress-common-lines -i -w -y -B <(cut -f 3,4,6,10 $DATADIR/$FILE.headless) $DATADIR/$FILE.recovered | wc -l
# 	# only compare chr, offs, cigar strings
# 	diff --suppress-common-lines -i -w -y -B <(cut -f 3,4,6 $DATADIR/$FILE) $DATADIR/$FILE.recovered | wc -l
# fi