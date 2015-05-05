#ifndef MERGED_EDITS_BUFFER_H
#define MERGED_EDITS_BUFFER_H

#include <cassert>
#include <memory>

#include "InputBuffer.hpp"
#include "RefereeUtils.hpp"

class MergedEditsStream {
public:
	MergedEditsStream(string const & file_name, uint8_t & read_len) {
		edits_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".edits") );
		read_len = (uint8_t)edits_in->getNextByte();
	}

	~MergedEditsStream() {
		// cerr << "Observed " << alignment_index << " alignments, of them " << have_edits << " had edits" << endl;
		// assert(alignment_index == alignments_expected);
	}

	vector<uint8_t> getEdits() {
		auto e = edits_in->getNextNBytes(num_edit_bytes);
		assert(e.size() == num_edit_bytes);
		return e;
	}

	bool hasEdits() {
		return num_edit_bytes > 0;
	}

	int next() {
		// if (num_edit_bytes) {
			// edits_in->popNBytes(num_edit_bytes);
		// }
		alignment_index++;
		if ( !edits_in->hasMoreBytes() ) {
			return END_OF_STREAM;
		}
		num_edit_bytes = (uint8_t)edits_in->getNextByte();
		have_edits += (num_edit_bytes > 0);
		return alignment_index;
	}

private:
	shared_ptr<InputBuffer> edits_in;

	uint8_t num_edit_bytes = 0;

	size_t alignment_index = 0;
	
	size_t have_edits = 0;

	// int alignments_expected = 0;


	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	inline shared_ptr<queue<string>> parseEdits(string const & fname, 
		char & read_len) {
		shared_ptr<queue<string>> edits(new queue<string>());
		string line;
		ifstream f_in(fname, ios::binary | ios::in);
		if (!f_in) {
			cerr << "[ERROR] Required file not found: " << fname << endl;
			exit(0);
		}

		auto fileSize = getFileLength(f_in);
		auto bytes_read = 0;
		char rlen = 0;
		f_in.read(reinterpret_cast<char *>(&rlen), 1);
		cerr << "Uniform read length: " << (int)rlen << endl;
		read_len = rlen;
		bytes_read++;
		unsigned char num_bytes = 0;

		while (bytes_read < fileSize) {
			// read one byte containing a number of edits to expect after it
			f_in.read(reinterpret_cast<char *>(&num_bytes), 1);
			// cerr << (int)num_bytes << "/'" << num_bytes << "' ";
			bytes_read++;
			assert(num_bytes > 0);
			string edit_str;
			edit_str.resize(num_bytes);
			f_in.read(&edit_str[0], num_bytes);
			edits->push(edit_str);
			bytes_read += num_bytes;
		}
		f_in.close();
		cerr << "Bytes read: " << bytes_read << endl;
		cerr << "Edit strings: " << edits->size() << endl;

		return edits;
	}
};

#endif
