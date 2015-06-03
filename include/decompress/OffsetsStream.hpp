#ifndef OFFSETS_BUFFER_H
#define OFFSETS_BUFFER_H

#include <memory>
#include "decompress/InputBuffer.hpp"


class OffsetsStream {
private:

	shared_ptr<InputBuffer> offsets_in;
	int current_multiplier = 0;
	int current_offset = 0;
	int delta = 0;
	int current_transcript = -1;

public:
	// OffsetsStream(string & file_name) {
	// 	offsets_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".offs") );
	// }

	OffsetsStream(shared_ptr<InputBuffer> ib): offsets_in(ib) { }

	////////////////////////////////////////////////////////////////////////
	pair<int, unsigned long> seekToBlockStart(int ref_id, int const start_coord, int const end_coord) {
		bool is_transcript_start = false;
		pair<int,unsigned long> p = offsets_in->loadOverlappingBlock(ref_id, start_coord, end_coord, is_transcript_start);
		auto block_start_coord = p.first;
		unsigned long block_num_align = p.second;
		// cerr << "2: " << block_num_align << endl;
		if (block_start_coord < 0) {
			cerr << "[ERROR] Could not navigate to the begining of the interval" << endl;
			exit(1);
		}
		current_multiplier = 0;
		current_offset = block_start_coord;
		// cerr << "PIZZA " << current_offset << endl;
		delta = 0;	// is delta set correctly?
		if (is_transcript_start) {
			// consumes bytes describing the transcript ID
			int ref_id = getNextTranscript();
			cerr << "current transcript: " << ref_id << endl;
		}
		// if actual_start_coord >= start_coord -- we did not have any alignment before this
		// and should start with actual_start_coord
		// if actual_start_coord < start_coord -- should handle this case in Decompressor and 
		// seek to the desired coordinate along w/ other data streams (clips, quals, edits, etc)
		p.first = getCurrentOffset();
		return p;
	}

	////////////////////////////////////////////////////////////////////////
	bool hasMoreOffsets() {
		if (!offsets_in->hasMoreBytes() ) return false;
		return true;
	}

	////////////////////////////////////////////////////////////////////////
	int getNextTranscript() {
		if ( !offsets_in->hasMoreBytes() ) return END_OF_STREAM;
		current_offset = 0;

		vector<uint8_t> chunk;
		char c = offsets_in->getNextByte();
		while (c != ' ' && offsets_in->hasMoreBytes()) {
			chunk.push_back(c);
			c = offsets_in->getNextByte();
		}
		if (chunk.size() == 0 || c==0) {
			cerr << "getting next transcript: no data" << endl;
			current_transcript = -1;
			return END_OF_STREAM;
		}
		current_transcript = stoi( string(chunk.begin(), chunk.end() ) );
		return current_transcript;
	}

	////////////////////////////////////////////////////////////////////////
	int getCurrentOffset() {
		return current_offset;
	}

	////////////////////////////////////////////////////////////////////////
	int getCurrentTranscript() {
		return current_transcript;
	}

	////////////////////////////////////////////////////////////////////////
	//
	// Returns a 0-based offset
	// Expected format:
	// <trans_id> <offset>:<multiplier> <offset> \n
	// <trans_id> <offset> \n
	//
	////////////////////////////////////////////////////////////////////////
	int getNextOffset() {
		if (current_transcript < 0) {
			cerr << "transcript < 0 " << current_transcript << endl;
			return END_OF_TRANS;
		}
		else if (current_multiplier > 0) {
			current_multiplier--; // used up one of the copies of this read
			return current_offset;
		}
		else {
			if ( !offsets_in->hasMoreBytes() ) return END_OF_STREAM;
			vector<uint8_t> chunk;
			char c = offsets_in->getNextByte();
			if (c == 0) {
				return END_OF_STREAM;
			}
			if (c == '\n') { // that's it for this transcript, moving on to the next
				cerr << "end of line" << endl;
				current_transcript = -1;
				return END_OF_TRANS;
			}
			while (c != ' ' && c != ':' && c != '\n' && offsets_in->hasMoreBytes() ) {
				chunk.push_back( c );
				c = offsets_in->getNextByte();
			}
			if (c == '\n') {
				cerr << "end of line" << endl;
				current_transcript = -1;
			}
			// delta from the previous absolute offset
			delta = stoi( string(chunk.begin(), chunk.end() ) );
			current_offset += delta;
			// now parse the multiplier if it exists
			if (c == ':') {
				chunk.clear();
				c = offsets_in->getNextByte();
				while (c != ' ' && c != '\n' && offsets_in->hasMoreBytes()) {
					chunk.push_back(c);
					c = offsets_in->getNextByte();
				}
				if (c == '\n') {
					cerr << "end of line" << endl;
					current_transcript = -1;
				}
				current_multiplier = stoi( string(chunk.begin(), chunk.end() ) );
				current_multiplier--;	// will return this offset once right now
			}
			return current_offset;
		}
	}
};

#endif