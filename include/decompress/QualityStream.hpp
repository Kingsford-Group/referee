#ifndef QUAL_STREAM_HPP
#define QUAL_STREAM_HPP


#include <memory>
#include "decompress/InputBuffer.hpp"


class QualityStream : public InputStream {
	// TODO: set up other buffers

public:

	// TODO: add buffers for clusters
	QualityStream(shared_ptr<InputBuffer> memb) : InputStream(memb) {
		// TODO: set up other buffers
	}

	// char *
	string getNextQualVector() {
		// TODO
		return "***";
	}


};

#endif