#ifndef OFFSETS_BUFFER_H
#define OFFSETS_BUFFER_H

#include <memory>
#include "decompress/InputBuffer.hpp"


class OffsetsStream {
public:
	// OffsetsStream(string & file_name) {
	// 	offsets_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".offs") );
	// }

	OffsetsStream(shared_ptr<InputBuffer> ib): offsets_in(ib) { }

	void seekTo(int ref_id, int start_coord) {
		bool is_transcript_start = false;
		auto actual_start_coordinate = offsets_in->loadOverlappingBlock(ref_id, start_coord, is_transcript_start);
		if (actual_start_coordinate < 0) {
			cerr << "[ERROR] Could not navigate to the begining of the interval" << endl;
			exit(1);
		}
		seekToCoord(ref_id, actual_start_coordinate, start_coord, is_transcript_start);
	}

	bool hasMoreOffsets() {
		if (!offsets_in->hasMoreBytes() ) return false;
		return true;
	}

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
			current_transcript = -1;
			return END_OF_STREAM;
		}
		current_transcript = stoi( string(chunk.begin(), chunk.end() ) );
		return current_transcript;
	}

	int getCurrentOffset() {
		return current_offset;
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
			return END_OF_TRANS;
		}
		else if (current_multiplier > 0) {
			current_multiplier--; // used up one of the copies of this read
			// current_offset += delta;
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
				current_transcript = -1;
				return END_OF_TRANS;
			}
			while (c != ' ' && c != ':' && c != '\n' && offsets_in->hasMoreBytes() ) {
				chunk.push_back( c );
				c = offsets_in->getNextByte();
			}
			if (c == '\n') {
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
					current_transcript = -1;
				}
				current_multiplier = stoi( string(chunk.begin(), chunk.end() ) );
				current_multiplier--;	// will return this offset once right now
			}
			return current_offset;
		}
	}

private:
	shared_ptr<InputBuffer> offsets_in;
	int current_multiplier = 0;
	int current_offset = 0;
	int delta = 0;
	int current_transcript = -1;


	void seekToCoord(int ref_id, int actual_start_coord, int start_coord, 
		bool const is_transcript_start) {
		cerr << "Seek to coord " << actual_start_coord << " " << start_coord << endl;
		current_multiplier = 0;
		delta = 0;
		if (is_transcript_start) {
			// consumes bytes describing the transcript ID
			getNextTranscript();
		}
		// if actual_start_coord >= start_coord -- we did not have any alignment before this
		// and should start with actual_start_coord
		// TODO
		// if actual_start_coord < start_coord -- should handle this case in Decompressor and 
		// seek to the desired coordinate along w/ other data streams (clips, quals, edits, etc)
	}
};

#endif