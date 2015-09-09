## Referee -- Rapid Separable Compression for Sequence Alignments

##### Dependencies

[Staden io_Lib](http://sourceforge.net/projects/staden/files/io_lib/)

[LZLib](http://www.nongnu.org/lzip/lzlib.html)

##### Installation



##### Usage

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
