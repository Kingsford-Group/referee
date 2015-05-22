#ifndef QUAL_STREAM_HPP
#define QUAL_STREAM_HPP


#include <memory>
#include "decompress/InputBuffer.hpp"


class QualityStream {
public:
	// QualityStream(string & file_name, int K_c) {
	// 	// open a stream for cluster membership
	// 	membership_stream = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".k=" + to_string(K_c) + ".membership"));
	// 	// TODO: open multiple streams for clusters and prefices/suffices
	// 	// quals_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".quals") );
	// }

	// TODO: add buffers for clusters
	QualityStream(shared_ptr<InputBuffer> memb) : membership_stream(memb) {}

	// char *
	string getNextQualVector() {
		// TODO
		return "***";
	}

private:
	shared_ptr<InputBuffer> membership_stream;
	shared_ptr<InputBuffer> quals_in;
};

#endif