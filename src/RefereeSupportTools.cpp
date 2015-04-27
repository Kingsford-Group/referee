/* Some useful tools: depth of coverage, number of edits */
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <cassert>

#include "decompress/OffsetsStream.hpp"

using namespace std;

////////////////////////////////////////////////
// does not account for splicing events
// are clipped regions part of the coverage?
////////////////////////////////////////////////
double depth(string & fname, size_t & genome_len, int const read_len) {
	// initialize vector long enough for a bacterial genome
	vector<uint32_t> covered_bases(6 * pow(10, 6), 0);

	cerr << "Max size: " << covered_bases.max_size() << endl;
	OffsetsStream offs(fname);

	// TODO: fill out covered_bases
	int ref_id = offs.getNextTranscript();
	uint32_t trans_offset = 0, last_offset = 0, huh = 0;

	while ( offs.hasMoreOffsets() ) {
		uint32_t offset = offs.getNextOffset();
		// switching to the next reference string
		if (offset == END_OF_TRANS) {
			trans_offset += last_offset + read_len;
			ref_id = offs.getNextTranscript();
			if (ref_id == END_OF_STREAM) {
				cerr << "Done" << endl;
				break;
			}
			cerr << "ref=" << ref_id << " ";
		}
		// done with everything
		else if (offset == END_OF_STREAM) {
			huh++;
		}
		// process the read
		else {
			if (covered_bases.size() < trans_offset + offset + read_len) {
				// add a million bases
				// cerr << covered_bases.capacity() << endl;
				// covered_bases.resize(trans_offset + offset + read_len + pow(10, 6), 0);
				uint32_t s = covered_bases.size();
				for (int i = 0; i < trans_offset + offset + read_len + pow(10,6) - s; i++)
					covered_bases.push_back(0);
				// cerr <<covered_bases.size() << " vs " << (trans_offset + offset + read_len ) << "| ";
				assert(covered_bases.size() >= trans_offset + offset + read_len );
			}
			// add 1 to all bases that this read covers
			for (int i = 0; i < read_len; i++) {
				assert(covered_bases.size() > trans_offset + offset + i);
				covered_bases[trans_offset + offset + i]++;
			}
			// remember the last offset
			last_offset = offset;
		}
	}

	cerr << "huh? " << huh << endl;

	uint32_t sum = accumulate(covered_bases.begin(), covered_bases.end(), (uint32_t)0);
	uint32_t nz = accumulate(covered_bases.begin(), covered_bases.end(), (uint32_t)0, 
		[](uint32_t & partial_sum, uint32_t & base){
			if (base > 0) {
				partial_sum++;
			}
			return partial_sum;
	});
	genome_len = nz;
	// cerr << sum << ", " << nz << endl;
	if (nz > 0) return sum / (double)nz;
	else return 0;
}

////////////////////////////////////////////////
int total_edits(string & fname) {
	string full_name = fname + ".edits.lz";
	ifstream edits_in(full_name);
	if (!edits_in) {
		cerr << "[ERROR] Could not open file. Terminating. " << full_name << endl;
		exit(1);
	}
	int edit_cnt = 0;
	edits_in.close();
	return edit_cnt;
}


////////////////////////////////////////////////
//
////////////////////////////////////////////////
int main(int argc, char * argv []) {
	if (argc < 4) {
        cerr << "Not enough arguments" << endl;
        exit(1);
    }
    string task = argv[1];
    string fname = argv[2];
    int read_len = stoi(argv[3]);

    if (task.compare("depth") == 0) {
    	// compute depth of coverage normalizing by # of covered bases
    	size_t g = 0;
    	double d = depth(fname, g, read_len);
    	cerr << "Average depth of coverage over covered bases: " << d << 
    		" (covered bases: " << g << ")" << endl;
    }
    else if (task.compare("edits") == 0) {
    	// compute the total number of edits (clips, mm, indels, splices)
    	total_edits(fname);
    }
}