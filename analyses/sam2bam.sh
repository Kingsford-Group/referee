#FILE=ERR005143_1.bam
#time samtools view -h $FILE -o $FILE.sam
DIR=/mnt/scratch0/dfilippo/aligned
# for FILE in SRR1294122 SRR445718 # MiSeq_Ecoli_DH10B_110721_PF # SRR1294122 # SRR445718 # K562_cytosol_LID8465_TopHat_v2 NA12878_S1
for FILE in NA12878_S1 # SRR445718 SRR1294122 # NA12878_S1 # MiSeq_Ecoli_DH10B_110721_PF # K562_cytosol_LID8465_TopHat_v2 NA12878_S1
do
        echo '------------'
        echo $FILE
        # sam to bam
        if [ ! -e $DIR/$FILE.bam ] 
        then
        	echo "convert from sam to unsorted bam"
        	time samtools view -b -o $DIR/$FILE.bam $DIR/$FILE.sam
        fi
        echo "sort the bam"
        time samtools sort -l 9 -@ 7 $DIR/$FILE.bam $DIR/$FILE.sorted
        # rm -rf $DIR/$FILE.bam

        echo "create index"
        time samtools index -b $DIR/$FILE.sorted.bam
        echo "convert from sorted bam to sorted sam"
        time samtools view -h $DIR/$FILE.sorted.bam > $DIR/$FILE.sorted.sam
		# rm $DIR/$FILE.sorted.bam
        ls -l $DIR/$FILE*
done

echo "Cha-cha!"
echo "If all is correct, remove original sam"
