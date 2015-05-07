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
DIR=/data/referee/aligned
FILE=SRR445718.2m.sam
# FILE=SRR445718.sam
#FILE=SRR1294122.sam
#DIR=.
# FILE=test.sam
# FILE=P_aeruginosa_PAO1
# PAERUG_GENOME=/data/genomes/bacterial/Pseudomonas_aeruginosa_PAO1_uid57945/NC_002516.fna
# HUMAN_GENOME=/mnt/scratch0/dfilippo/genomes/deez/all_chromosomes_hg_19.fa
# GENOME=$DIR/human_genome.fa
GENOME=$DIR/human_chr10_37.2.fa
#FILE=K562_cytosol_LID8465_TopHat_v2.sam
# FILE=NA12878_S1.sorted.sam
# rm -f $DIR/$FILE.offs.lz
# rm -f $DIR/$FILE.edits.lz
# rm -f $DIR/$FILE.left_clip.lz
# rm -f $DIR/$FILE.right_clip.lz
# rm -f $DIR/$FILE.flags.lz
# rm -f $DIR/$FILE.ids.lz
# rm -f $DIR/$FILE.quals.lz
# rm -f $DIR/$FILE.opt.lz
# rm -f $DIR/$FILE.unaligned.lz
# rm -f $DIR/$FILE.k\=*

# HUMAN_GENOME=$DIR/human_genome.fa
# HUMAN_GENOME=$DIR/../genomes/deez/all_chromosomes_hg_19.fa

# TIMEOPT="-v"
TIMEOPT="-lp"

/usr/bin/time $TIMEOPT $BIN -c $DIR/$FILE -t 10 -r $GENOME # > $DIR/$FILE.seq_comp.log 2>&1
# ls -l $DIR/$FILE.*
# time plzip -vf $DIR/$FILE.k\=* 2> $DIR/$FILE.plzip
# python python/parse_plzip_output.py $DIR/$FILE.plzip
# $BIN -c $DIR/$FILE 10
# ls -l $DIR/$FILE.*
# /usr/bin/time -lp $BIN -c $DIR/P_aeruginosa_PAO1.10mil.sam

# DECOMPRESS=1
# if [ $DECOMPRESS -eq 1 ]
# then
# 	plzip -fd $DIR/$FILE.offs.lz 2> /dev/null
# 	plzip -fd $DIR/$FILE.edits.lz 2> /dev/null
# 	plzip -fd $DIR/$FILE.*clip.lz 2> /dev/null
# 	plzip -fd $DIR/$FILE.has_edits.lz 2> /dev/null
# 	plzip -fd $DIR/$FILE.*.lz 2> /dev/null
# 	# /usr/bin/time $TIMEOPT $BIN -d $DIR/$FILE -r $GENOME # > $DIR/$FILE.decomp.log 2>&1
# 	# only compare chr, offs, cigar strings
# 	# diff --suppress-common-lines -i -w -y -B <(cut -f 3,4 $DIR/$FILE) $DIR/$FILE.recovered | wc -l

# 	# compare IDs
# 	# echo "Compare read IDs (only for the aligned reads)"
# 	# diff --suppress-common-lines -i -w -y -B <(gawk -F't' '{ if (and($2,4) == 0) print $1 }' SRR445718.2m.sam.headless) <(gawk -F't' '{ if (and($2,4) == 0) print $1 }' SRR445718.2m.sam.recovered) | wc -l
# 	# achieved 0 differences

# 	echo "Compare first 10 columns"
# 	FIN=$DIR/$FILE.headless
# 	FOUT=$DIR/$FILE.recovered
# 	PARAMS="--suppress-common-lines -i -w -y -B"
# 	COLUMNS="\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10"
# 	diff $PARAMS <(gawk -F'\t' '{ if (and($2,4) == 0) print $COLUMNS }' $FIN) <(gawk -F'\t' '{ if (and($2,4) == 0) print $COLUMNS }' $FOUT) | wc -l
# fi