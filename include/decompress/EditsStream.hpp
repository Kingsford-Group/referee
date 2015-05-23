#ifndef EDITS_BUFFER_H
#define EDITS_BUFFER_H

#include <cassert>
#include <memory>

#include "InputBuffer.hpp"

class EditsStream {
public:
	// EditsStream(string const & file_name) {
	// 	edits_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".edits") );
	// 	if (!edits_in->opened()) {
	// 		cerr << "[ERROR] Required file missing. Terminating." << endl;
	// 		exit(1);
	// 	}
	// 	// first 2 bytes contain the read length
	// 	read_len = edits_in->getNextByte();
	// 	bytes_read++;

	// 	has_edits_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".has_edits") );
	// 	if (!has_edits_in->opened()) {
	// 		cerr << "[ERROR] Required file missing. Terminating." << endl;
	// 		exit(1);
	// 	}
	// 	// first 4 bytes hold the number of alignments total
	// 	uint8_t c = has_edits_in->getNextByte();
	// 	alignments_expected |= ( c << 24) ;
	// 	c = has_edits_in->getNextByte();
	// 	alignments_expected |= ( c << 16) ;
	// 	c = has_edits_in->getNextByte();
	// 	alignments_expected |= ( c << 8) ;
	// 	c = has_edits_in->getNextByte();
	// 	alignments_expected |= ( c ) ;

	// 	has_edit_byte = has_edits_in->getNextByte();
	// 	// cerr << "First has edit byte: " << (int)has_edit_byte << " i=" << (int)i << endl;
	// }


	EditsStream(shared_ptr<InputBuffer> e, shared_ptr<InputBuffer> h): 
		edits_in(e), 
		has_edits_in(h) {}

	~EditsStream() {
		cerr << "Expected: " << alignments_expected << " observed: " << alignment_count << endl;
		// assert(alignment_count == alignments_expected);
	}

	short getReadLen() { return read_len; }

	size_t getAlignmentCount() {return alignments_expected; }

	vector<uint8_t> getEdits() {
		// cerr << "Bytes read before polling more: " << bytes_read << " vs. " << edits_in->current_index << endl;
		uint8_t num_edit_bytes = edits_in->getNextByte();
		bytes_read++;
		auto e = edits_in->getNextNBytes(num_edit_bytes);
		bytes_read += num_edit_bytes;
		// cerr << "Edit seq len: " << (int)num_edit_bytes << " (" << bytes_read << " vs. " << edits_in->current_index  << ")" << endl;
		assert(e.size() == num_edit_bytes);
		return e;
	}


	bool hasEdits() {
		// cerr << "Edit byte: " << (int)has_edit_byte << " i=" << (int)i << endl;
		return (has_edit_byte >> i) & 1;
	}

	int next() {
		alignment_count++;
		i--;
		// cerr << "NEXT i=" << (int)i << endl;
		if (i < 0) {
			if ( !has_edits_in->hasMoreBytes() ) {
				return END_OF_STREAM;
			}
			has_edit_byte = has_edits_in->getNextByte();
			i = sizeof(has_edit_byte) * 8 - 1;
		}
		return SUCCESS;
	}

private:
	shared_ptr<InputBuffer> edits_in;

	int bytes_read = 0;

	shared_ptr<InputBuffer> has_edits_in;

	size_t alignment_count = 0;

	short read_len = 0;

	char i = sizeof(has_edit_byte) * 8;

	uint8_t has_edit_byte;

	size_t alignments_expected = 0;
};

#endif