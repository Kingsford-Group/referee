#ifndef FLAGS_STREAM_HPP
#define FLAGS_STREAM_HPP


#include <memory>
#include "decompress/InputBuffer.hpp"


class FlagsStream {
public:
	// FlagsStream(string & file_name/*, TranscriptsStream const & transcripts*/) {
	// 	flags_in = shared_ptr<InputBuffer>(new InputBuffer(file_name + ".flags") );
	// 	// blah = transcripts.getFlagDictionary();
	// }

	FlagsStream(shared_ptr<InputBuffer> in) : flags_in(in) {}

	vector<int> getNextFlagSet() {
		vector<int> loc_flags;

		if ( !flags_in->hasMoreBytes() ) return loc_flags;

		string chunk;
		char c = flags_in->getNextByte();
		while (c != '\n' && flags_in->hasMoreBytes()) {
			if (c == ' ') {
				int value = stoi(chunk);
				// if mapped: dictionary[mapped_value]
				loc_flags.push_back(value);
				chunk = "";
			}
			else
				chunk.push_back(c);
			c = flags_in->getNextByte();
		}
		if (chunk.size() == 0 || c == 0) {
			return loc_flags;
		}
		return loc_flags;
	}

private:
	shared_ptr<InputBuffer> flags_in;
};

#endif
