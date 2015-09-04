## Referee -- Rapid Separable Compression for Sequence Alignments

##### Dependencies

Staden IO_Lib

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
