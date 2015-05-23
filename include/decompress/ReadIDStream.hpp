#ifndef READ_ID_STREAM_HPP
#define READ_ID_STREAM_HPP


#include <memory>
#include "decompress/InputBuffer.hpp"


class ReadIDStream {
public:
	// ReadIDStream(string & file_name) {
	// 	read_ids_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".ids") );
	// }

	ReadIDStream(shared_ptr<InputBuffer> in) : read_ids_in(in) {}

	int getNextID(string & id) {
		// exhausted the input stream
		if ( !read_ids_in->hasMoreBytes() ) return END_OF_STREAM;

		string chunk;
		char c = read_ids_in->getNextByte();
		while (c != '\n' && read_ids_in->hasMoreBytes()) {
			chunk.push_back(c);
			c = read_ids_in->getNextByte();
		}
		id = chunk;
		if (chunk.size() == 0 || c == 0) {
			return END_OF_STREAM;
		}
		return SUCCESS;
	}

private:
	shared_ptr<InputBuffer> read_ids_in;
};

#endif