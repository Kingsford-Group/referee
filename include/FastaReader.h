#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "unistd.h"
#include "fcntl.h"

#ifndef KSEQ
// KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(int, read)
#endif

class FastaReader {
    private:
    int fp; // file handler
    kseq_t *seq;
    int l;

    public:

    FastaReader(char * fname) {
        fp = open(fname, O_RDONLY); // TODO: check params
        seq = kseq_init(fp);
    }

    ~FastaReader() {
        kseq_destroy(seq);
	close(fp);
    }

    kseq_t * nextSequence() {
        l = kseq_read(seq);
        if (l < 0) return NULL;
        else return seq;
    }

};

#endif /* FASTA_READER_H */
