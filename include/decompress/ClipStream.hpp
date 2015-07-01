#ifndef CLIP_STREAM_H
#define CLIP_STREAM_H

#include <memory>

#include "decompress/InputStream.hpp"


class ClipStream : public InputStream {

	string current_clip;

public:

	ClipStream(shared_ptr<InputBuffer> ib): InputStream(ib) {}

	string peekNext() {
		// cerr << "current_clip: " << current_clip;
		return current_clip;
	}

	int getNext(string & clip) {
		if ( !data_in->hasMoreBytes() ) return END_OF_STREAM;
		string chunk;
		char c = data_in->getNextByte();
		while (c != '\n' && data_in->hasMoreBytes()) {
			chunk.push_back(c);
			c = data_in->getNextByte();
		}
		clip = chunk;
		// cerr << "Clip: " << clip << endl;
		current_clip = chunk;
		if (chunk.size() == 0 || c == 0) {
			return END_OF_STREAM;
		}
		return SUCCESS;
	}

};

#endif