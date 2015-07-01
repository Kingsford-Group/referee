#ifndef READ_ID_STREAM_HPP
#define READ_ID_STREAM_HPP


#include <memory>
#include "decompress/InputStream.hpp"


class ReadIDStream : public InputStream {

public:

	ReadIDStream(shared_ptr<InputBuffer> buf) : InputStream(buf) {}

	string getNextID(int & status) {
		// exhausted the input stream
		// cerr << "IDs: " << data_in << endl;
		if ( !data_in->hasMoreBytes() ) {
			status = END_OF_STREAM;
			return "";
		}

		string chunk = "";
		char c = data_in->getNextByte();
		
		while (c != '\n' && data_in->hasMoreBytes()) {
			// cerr << c;
			chunk.push_back(c);
			c = data_in->getNextByte();
		}
		// cerr << "cha" << endl;
		if (chunk.size() == 0 || c == 0) {
			status = END_OF_STREAM;
			return "";
		}
		status = SUCCESS;
		return chunk;
	}


};

#endif