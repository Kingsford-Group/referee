#ifndef FLAGS_STREAM_HPP
#define FLAGS_STREAM_HPP


#include <memory>
#include "decompress/InputBuffer.hpp"


class FlagsStream {
public:
	FlagsStream(string & file_name/*, TranscriptsStream const & transcripts*/) {
		flags_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".flags") );
		// blah = transcripts.getFlagDictionary();
	}

	int getNextFlagSet(vector<int> & flags) {
		// TODO
		vector<int> flags;

		if ( !flags_in->hasMoreBytes() ) return END_OF_STREAM;

		string chunk;
		char c = flags_in->getNextByte();
		while (c != '\n' && flags_in->hasMoreBytes()) {
			if (c == ' ') {
				int value = stoi(chunk);
				// if mapped: dictionary[mapped_value]
				flags.push_back(value);
				chunk = "";
			}
			else
				chunk.push_back(c);
			c = flags_in->getNextByte();
		}
		id = chunk;
		if (chunk.size() == 0 || c == 0) {
			return END_OF_STREAM;
		}
		return SUCCESS;
	}

private:
	shared_ptr<InputBuffer> flags_in;
};

#endif