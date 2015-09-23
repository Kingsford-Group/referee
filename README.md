## Referee: Rapid, Separable Compression for Sequence Alignments

Referee is a command-line tool that takes sequence alignment SAM files and compresses them in a lossless manner. Referee achieves significant improvements in sequence compression: it can compress 187mln bases down to 135Mb of space, an equivalent of 0.06 bits required to encode a single base. Its streaming clustering algorithm for quality values improves over state of the art tools while allowing random access to the compressed data. Additionally, you can operate on the compressed streams without having to decompress your data. For example, calculating sequencing depth on 36mln alignments of human embryonic stem cells (SRR445718) takes about 8s (for comparison, samtools takes 19m on the same dataset). Finally, you can decompress data streams independently and skip on the downloads for read IDs, quality values, optional SAM flags, and so on. For details, see the corresponding publication below.


#### Dependencies

[Staden io_Lib](http://sourceforge.net/projects/staden/files/io_lib/)

[LZLib](http://www.nongnu.org/lzip/lzlib.html)

#### Installation



#### Usage

To compress:

```
	referee [options] -r reference.fa alignments.sam
```

To decompress:

```
	referee -d [options] -r reference.fa alignments.sam
```

Options:

	-t=N                 number of threads

	--seqOnly            encode sequencing data only

	--discardSecondary   discard secondary alignments

	view chrK:L-M        retrieve data from interval [L,M) on chromosome K

	-h, --help           this help


#### Cite

Filippova, D., Kingsford, C. Rapid, Separable Compression Enables Fast Analyses of Sequence Alignments. ACM BCB, Atlanta. 2015.
