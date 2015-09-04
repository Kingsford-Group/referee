#ifndef INPUT_STREAM_LIB_H
#define INPUT_STREAM_LIB_H

#include "decompress/InputBuffer.hpp"

class InputStream {
protected:

	shared_ptr<InputBuffer> data_in;

public:

	InputStream(shared_ptr<InputBuffer> in) : data_in(in) {}

	pair<int, unsigned long> seekToBlockStart(int const ref_id, int const start_coord, int const end_coord) {
		bool t = false;
		auto start = data_in->loadOverlappingBlock(ref_id, start_coord, end_coord, t);
		if (start.first < 0) {
			cerr << "[ERROR] Could not navigate to the begining of the interval" << endl;
			exit(1);
		}
		// cerr << "loaded overlapping block" << endl;
		return start;
	}
};


#endif INPUT_STREAM_LIB_H