#ifndef EDITS_BUFFER_H
#define EDITS_BUFFER_H

#include <cassert>
#include <memory>

#include "InputBuffer.hpp"

class EditsStream {
private:
	shared_ptr<InputBuffer> edits_in;

	int bytes_read = 0;

	shared_ptr<InputBuffer> has_edits_in;

	size_t alignment_count = 0;

	short read_len = 0;

	char i = sizeof(has_edit_byte) * 8;

	uint8_t has_edit_byte;

	size_t alignments_expected = 0;

	////////////////////////////////////////////////////////////////////////////
	// 
	////////////////////////////////////////////////////////////////////////////
	pair<int,unsigned long> syncEditStreams(pair<int, unsigned long> & edits_start, 
		pair<int, unsigned long> & has_edits_start) {
		int edit_start_coord = edits_start.first;
		unsigned long edit_num_al = edits_start.second;
		// int has_edit_start_coord = has_edits_start.first;
		unsigned long has_edits_num_al = has_edits_start.second;

		// should never happen since we request that has_edits genomic coordinate is
		// less than or equal to the earliest edits genomic coordinate
		// while (edit_num_al < has_edit_num_al) {
		// }
		cerr << edit_num_al << " vs " << has_edits_num_al << endl;
		assert(edit_num_al >= has_edits_num_al);

		while (has_edits_num_al < edit_num_al) {
			// need off to seek forward to edit_start
			has_edits_in->getNextByte();
			has_edits_num_al++;
		}
		cerr << "has_edits_num_al: " << has_edits_num_al << endl;
		// can not seek to the target genomic coord since we do not know the mapping between 
		// genomic coords and edits/has_edits stream
		// only the offset stream knows this mapping
		return make_pair(edit_start_coord, has_edits_num_al);
	}

public:

	////////////////////////////////////////////////////////////////////////////
	EditsStream(shared_ptr<InputBuffer> e, shared_ptr<InputBuffer> h): 
		edits_in(e), 
		has_edits_in(h) {
		}

	////////////////////////////////////////////////////////////////////////////
	~EditsStream() {
		cerr << "Expected: " << alignments_expected << " observed: " << alignment_count << endl;
		// assert(alignment_count == alignments_expected);
	}

	////////////////////////////////////////////////////////////////////////////
	// Loads the first block of data that overlapping the coordinate and syncs the underlying
	// data streams to point at the same starting alignment at genomic coordinate less than or
	// equal to the input (target) coordinate
	////////////////////////////////////////////////////////////////////////////
	pair<int, unsigned long> seekToBlockStart(int const ref_id, int const start_coord, int const end_coord) {
		bool transcript_start = false;
		auto edits_start = edits_in->loadOverlappingBlock(ref_id, start_coord, end_coord, transcript_start);
		if (edits_start.first < 0) {
			cerr << "[ERROR] Could not navigate to the begining of the interval" << endl;
			exit(1);
		}
		// pass the start coordinate for the edits block -- this way we can still sync
		auto has_edits_start = has_edits_in->loadOverlappingBlock(ref_id, edits_start.first, end_coord, transcript_start);
		if (has_edits_start.first < 0) {
			cerr << "[ERROR] Could not navigate to the begining of the interval" << endl;
			exit(1);
		}
		// sync these streams
		auto synced_coord = syncEditStreams(edits_start, has_edits_start);
		return synced_coord;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	short getReadLen() { return read_len; }

	//////////////////////////////////////////////////////////////////////////////////////////////
	size_t getAlignmentCount() {return alignments_expected; }

	//////////////////////////////////////////////////////////////////////////////////////////////
	vector<uint8_t> getEdits() {
		uint8_t num_edit_bytes = edits_in->getNextByte();
		bytes_read++;
		auto e = edits_in->getNextNBytes(num_edit_bytes);
		bytes_read += num_edit_bytes;
		assert(e.size() == num_edit_bytes);
		return e;
	}

	bool hasEdits() {
		return (has_edit_byte != 0);
	}

	// byte version 
	int next() {
		alignment_count++;
		if ( !has_edits_in->hasMoreBytes() ) {
			return END_OF_STREAM;
		}
		has_edit_byte = has_edits_in->getNextByte();
		return SUCCESS;
	}


};

#endif