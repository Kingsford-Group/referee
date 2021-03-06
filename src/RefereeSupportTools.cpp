/* Some useful tools: depth of coverage, number of edits */
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <cassert>
#include <deque>

#include "decompress/OffsetsStream.hpp"
#include "decompress/EditsStream.hpp"

using namespace std;

////////////////////////////////////////////////
// does not account for splicing events
// are clipped regions part of the coverage?
// todo: can implement w/o keeping the whole array around
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

	size_t sum = accumulate(covered_bases.begin(), covered_bases.end(), (size_t)0);
	size_t nz = accumulate(covered_bases.begin(), covered_bases.end(), (size_t)0,
		[](size_t & partial_sum, uint32_t & base){
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
// does not account for splicing events
// are clipped regions part of the coverage?
// does not keep the whole array for covered bases around (low-mem and faster)
////////////////////////////////////////////////
double rolling_depth(string & fname, size_t & genome_len, int const read_len) {
	// initialize vector long enough for a bacterial genome
	deque<int> covered_bases(read_len, 0);
	size_t sum = 0;
	OffsetsStream offs(fname);
	int ref_id = offs.getNextTranscript();
	int last_offset = 0, huh = 0;
	cerr << "ref=" << ref_id+1 << " ";
	while ( offs.hasMoreOffsets() ) {
		int offset = offs.getNextOffset();
		if (last_offset == 0) {
			last_offset = offset;
		}
		// switching to the next reference string
		if (offset == END_OF_TRANS) {
			last_offset = 0;
			// add up counts in covered_bases; reset to 0's
			for (int i = 0; i < covered_bases.size(); i++) {
				sum += covered_bases[i];
				covered_bases[i] = 0;
				genome_len++;
			}
			ref_id = offs.getNextTranscript();
			if (ref_id == END_OF_STREAM) {
				cerr << "Done" << endl;
				break;
			}
			cerr << "ref=" << ref_id+1 << " ";
		}
		// done with everything
		else if (offset == END_OF_STREAM) {
			huh++;
		}
		// process the read
		else {
			// fold the items on the left that are out of range
			// cerr << "off=" << offset << " last_off=" << last_offset <<  " ";
			if (offset - last_offset > read_len) {
				last_offset = offset - read_len;
			}
			for (int i = 0; i < offset - last_offset; i++ ) {
				sum += covered_bases.front(); // fold
				covered_bases.pop_front();	// pop
				covered_bases.push_back(0);	// push a new count on the end
				genome_len++;
			}
			// add 1 to all bases that this read covers
			for (int i = 0; i < read_len; i++) {
				covered_bases[i]++;
			}
			// remember the last offset
			last_offset = offset;
		}
	}	
	if (genome_len > 0) return sum / (double)genome_len;
	else return 0;
}

////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////
size_t total_edits(string & fname) {
	uint8_t readlen;
	EditsStream edits(fname);
	size_t edit_cnt = 0;
	size_t total_alignments = edits.getAlignmentCount();
	short read_len = edits.getReadLen();
	// if counted inside the loop -- get 2x of actual alignment count
	size_t alignments_with_edits = 0;

	int return_value = 0;
	// size_t alignments = 0;
	while ( (return_value = edits.next() ) != END_OF_STREAM) { // advance to the next alignment)
		// alignments++;
		// total_alignments++;
		if (edits.hasEdits() ) {
			alignments_with_edits++;
			vector<uint8_t> edit_ops = edits.getEdits();
			for (int i = 0; i < edit_ops.size(); i++) {
				switch (edit_ops[i]) {
					case 'L':
					case 'R':
						// do nothing
					break;
					case 'E':
						i+=3;
					break;
					case 197:
						i += 4;
					break;
					case 'l':
					case 'r':
					case 'D':
					case 'V':
					case 'W':
					case 'X':
					case 'Y':
					case 'Z':
						i++;
					break;
				}
				edit_cnt++;
			}
		}
	}
	cerr << "Saw " << total_alignments << " total alignments; of them " << alignments_with_edits << " had edits." << endl;
	cerr << "Error rate: " << (long double)edit_cnt / (total_alignments * read_len) * 100 << "%" << endl;
	// TODO: total_alignments almost 2x that of actual number. figure it out
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
    	double d = rolling_depth(fname, g, read_len);
    	cerr << "Average depth of coverage over covered bases: " << d << 
    		" (covered bases: " << g << ")" << endl;
    }
    else if (task.compare("edits") == 0) {
    	// compute the total number of edits (clips, mm, indels, splices)
    	size_t edit_cnt = total_edits(fname);
    	cerr << "Total edits: " << edit_cnt << endl;
    }
}
