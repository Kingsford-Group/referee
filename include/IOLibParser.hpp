#ifndef IOLIB_PARSER_H
#define IOLIB_PARSER_H

#include <iostream>
#include <string>

extern "C" {
    #include "io_lib/scram.h"
    #include "io_lib/os.h"
    #undef max
    #undef min
}

using namespace std;

class IOLibParser {
public:
	IOLibParser (std::string const & file_name, int t) {
		// open file
	    fp_ = scram_open(file_name.c_str(), "r");
	    if (fp_ == NULL) {
	    	// failed to open the file
	    	cerr << "[ERROR] Can not open file: " << file_name << endl;
	    	cerr << "Please check that the file exists." << endl;
	    	failed_ = true;
	    	return;
	    }
	    scram_set_option(fp_, CRAM_OPT_NTHREADS, t /* threads */);

	    // read header
	    header_ = scram_get_header(fp_);
	    sam_hdr_incr_ref(header_);
	}

	bool failed() {return failed_;}

	bool read_next() {
		bool result = scram_get_seq(fp_, &the_read) >= 0;
		if (result) lines++;

		if (lines % 1000000 == 0) {
            std::cerr << lines / 1000000 << "mln ";
        }

		return result;
	}

	SAM_hdr* header() {
		return header_;
	}

	bam_seq_t* getRead() {
		return the_read;
	}

	void close() {
		std::cerr << "Read " << lines << " alignments" << std::endl;
		scram_close(fp_);
	}

private:
	int lines = 0;
	scram_fd* fp_ = NULL; 		// file pointer
	SAM_hdr* header_ = NULL;		// SAM file header
	bam_seq_t* the_read = NULL;	// pointer to an alignment
	bool failed_ = false;
};

#endif