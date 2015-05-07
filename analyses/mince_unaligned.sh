MINCE=/home/robp/mince/build/src
DIR=/mnt/scratch0/dfilippo/aligned
# make sure tbb is ON!
export LD_LIBRARY_PATH=/data/srameta/scripts-src/SalmonBeta-v0.2.2_ubuntu-14.04/lib/:$LD_LIBRARY_PATH
for file in SRR445718 # SRR1294122 P_aeruginosa_PAO1 K562_cytosol_LID8465_TopHat_v2
do
	FASTQ=$DIR/$file.sam.unaligned
	ls -l $FASTQ.lz
	plzip -dc -n 10 $FASTQ.lz > $FASTQ.fq
	/usr/bin/time -v $MINCE/mince -e -l U -r $FASTQ.fq -o $FASTQ.mince
done
