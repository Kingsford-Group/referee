#ifndef CLIP_STREAM_H
#define CLIP_STREAM_H

#include <memory>

#include "InputBuffer.hpp"


class ClipStream {
public:
	// ClipStream(string & file_name, char const * suffix) {
	// 	clips_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + suffix) );
	// }

	ClipStream(shared_ptr<InputBuffer> ib): clips_in(ib) {}

	pair<int, unsigned long> seekToBlockStart(int const ref_id, int const start_coord, int const end_coord) {
		cerr << "Clips buf: " << clips_in << endl;
		bool t = false;
		auto start = clips_in->loadOverlappingBlock(ref_id, start_coord, end_coord, t);
		cerr << "loaded overlapping block" << endl;
		if (start.first < 0) {
			cerr << "[ERROR] Could not navigate to the begining of the interval" << endl;
			exit(1);
		}
		return start;
	}


	string peekNext() {
		// cerr << "current_clip: " << current_clip;
		return current_clip;
	}

	int getNext(string & clip) {
		if ( !clips_in->hasMoreBytes() ) return END_OF_STREAM;
		string chunk;
		char c = clips_in->getNextByte();
		while (c != '\n' && clips_in->hasMoreBytes()) {
			chunk.push_back(c);
			c = clips_in->getNextByte();
		}
		clip = chunk;
		// cerr << "Clip: " << clip << endl;
		current_clip = chunk;
		if (chunk.size() == 0 || c == 0) {
			return END_OF_STREAM;
		}
		return SUCCESS;
	}

private:
	shared_ptr<InputBuffer> clips_in;
	string current_clip;
};

#endif