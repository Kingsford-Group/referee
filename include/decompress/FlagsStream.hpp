#ifndef FLAGS_STREAM_HPP
#define FLAGS_STREAM_HPP


#include <memory>
#include "decompress/InputBuffer.hpp"


class FlagsStream {

	shared_ptr<InputBuffer> flags_in;

	unordered_map<int,short> & r_flags_map;

	unordered_map<int,int> & r_mapq_map;

	unordered_map<int,int> & r_rnext_map;

	void translate(vector<int> & indexed_flags) {
		// flags
		if (r_flags_map.find(indexed_flags[0]) != r_flags_map.end())
			indexed_flags[0] = r_flags_map[indexed_flags[0]];
		// mapq
		if (r_mapq_map.find(indexed_flags[1]) != r_mapq_map.end())
			indexed_flags[1] = r_mapq_map[indexed_flags[1]];
		// rnext
		if (r_rnext_map.find(indexed_flags[2]) != r_rnext_map.end())
			indexed_flags[2] = r_rnext_map[indexed_flags[2]];
	}

public:

	// constructor
	FlagsStream(shared_ptr<InputBuffer> in, 
		unordered_map<int,short> & flags_map,
		unordered_map<int,int> & mapq_map,
		unordered_map<int,int> & rnext_map) : 
	flags_in(in), 
	r_flags_map(flags_map), 
	r_mapq_map(mapq_map), 
	r_rnext_map(rnext_map) {}

	// sync the stream to a specific coordinate
	// ref_id -- chromosome index
	// start_coord -- base pair address
	pair<int, unsigned long> seekToBlockStart(int const ref_id, int const start_coord, int const end_coord) {
		// cerr << "Flags buf: " << clips_in << endl;
		bool t = false;
		auto start = flags_in->loadOverlappingBlock(ref_id, start_coord, end_coord, t);
		cerr << "loaded overlapping block" << endl;
		if (start.first < 0) {
			cerr << "[ERROR] Could not navigate to the begining of the interval" << endl;
			exit(1);
		}
		return start;
	}

	short translateFlag(int value) {
		if (r_flags_map.find(value) != r_flags_map.end())
			value = r_flags_map[value];
	}

	// get next set of flags for an alignment
	vector<int> getNextFlagSet() {
		vector<int> loc_flags;
		if ( !flags_in->hasMoreBytes() ) return loc_flags;

		string chunk;
		char c = flags_in->getNextByte();
		while (c != '\n' && flags_in->hasMoreBytes()) {
			if (c == ' ') {
				int value = stoi(chunk);
				loc_flags.push_back(value);
				chunk = "";
			}
			else
				chunk.push_back(c);
			c = flags_in->getNextByte();
		}
		// add the last 'word' we observed before \n
		loc_flags.push_back(stoi(chunk));

		// translate values
		translate(loc_flags);

		if (chunk.size() == 0 || c == 0) {
			return loc_flags;
		}
		return loc_flags;
	}
};

#endif
