#ifndef QUAL_STREAM_HPP
#define QUAL_STREAM_HPP


#include <memory>
#include <numeric>
#include "decompress/InputBuffer.hpp"

#define END_OF_STREAM -2

class QualityStream : public InputStream {

	////////////////////////////////////////////////////////
	vector<int> seen_so_far;

	vector<shared_ptr<InputBuffer>> cores;

	vector<shared_ptr<InputBuffer>> prefixes;

	vector<shared_ptr<InputBuffer>> suffixes;

	shared_ptr<InputBuffer> other_qvs;

	// attempts to retrive the next line from the quality stream (core, prefix or suffix)
	string getLine(shared_ptr<InputBuffer> buf) {
		if ( !buf->hasMoreBytes() ) {
			cerr << "[ERROR] Requesting more data from quality stream, but stream is depleted." << endl;
			// exit(1);
			return "*";
		}
		string chunk;
		char c = buf->getNextByte();
		while (c != '\n') {
			chunk.push_back(c);
			if (!buf->hasMoreBytes()) break;
			c = buf->getNextByte();
		}
		return chunk;
	}

	int getNextMember() {
		if ( !data_in->hasMoreBytes() ) {
			cerr << "total algn seen: " << accumulate(seen_so_far.begin(), seen_so_far.end(), 0) << endl;
			for (int c : seen_so_far) cerr << c << " ";
				cerr << endl;
			cerr << "[ERROR] Requesting more data from (membership) quality stream, but stream is depleted." << endl;
			exit(1);
			return -1;
		}
		string chunk;
		char c = data_in->getNextByte();
		// chunk.push_back(c);
		while (c != ' ') {
			chunk.push_back(c);
			if (!data_in->hasMoreBytes()) break;
			c = data_in->getNextByte();
		}
		int value = stoi(chunk);
		return value;
	}

	/*
	Explode the compact string format to recover the original core

	Compacting code in compress/QualityCluster::writeCore(string const & core, GenomicCoordinate & currentCoord)

	# python pseudocode
	# for i in xrange(n):
	# 	symbol = line[i]
	# 	if symbol.isdigit():
	# 		# remap to a 130 range
	# 		symbol =  chr(130 + int(symbol) )
	# 	if symbol != 'J':
	# 		delta = i - prev_pos - 1
	# 		if delta > 0:
	# 			edit_str += str(delta)
	# 		edit_str += symbol
	# 		prev_pos = i
	*/
	string explodeString(const string & q_v) {
		// cerr << "explode string: " << q_v << endl;
		string core;
		bool had_low_quals = (q_v.find('!') != string::npos);
		if (had_low_quals) {
			// was not transformed
			return q_v;
		}
		else {
			// string was transformed
			for (int i = 0; i < q_v.size(); i++) {
				if ( (uint8_t)q_v[i] > 130) {
					// was a digit
					uint8_t d = q_v[i] - (uint8_t)130;
					// cerr << (char)d << " ";
					core.push_back(d);
				}
				else {
					string digits;
					while (isdigit(q_v[i]) ) {
						// cerr << q_v[i] << " ";
						digits.push_back(q_v[i]);
						i++;
					}
					int reps = 1;
					if (digits.size() > 0) {
						// cerr << "digits: " << digits << " ";
						reps = stoi(digits);
					}
					uint8_t symbol = q_v[i];
					// cerr << symbol << " ";
					for (int j = 0; j < reps; j++) {
						core.push_back(symbol);
					}
				}
			}
			// cerr << endl;
		}
		// cerr << "recovered: " << core << endl;
		return core;
	}

public:

	// add buffers for clusters
	QualityStream(shared_ptr<InputBuffer> memb, const string & path, 
		const unordered_map<string,shared_ptr<vector<TrueGenomicInterval>>> & all_intervals,
		const int buffer_size) : 
		InputStream(memb) {
			// cerr << "settting up quals stream" << endl;
		auto i = path.find(".membership");
		auto prefix = path.substr(0, i);
		// cerr << prefix << endl;
		bool exists = true;
		i = 1;
		while (exists) {
			string core_suf = ".quals." + to_string(i) + ".lz";
			string prefix_suf = ".quals." + to_string(i) + ".prefix.lz";
			string suffix_suf = ".quals." + to_string(i) + ".suffix.lz";
			
			if ( all_intervals.find(core_suf) == all_intervals.end() ) break;

			shared_ptr<InputBuffer> cluster_core(new InputBuffer(prefix + core_suf,
				all_intervals.find(core_suf)->second, buffer_size, 0) );
			// check that it opened
			// else exists = false;
			cores.push_back(cluster_core);
			shared_ptr<InputBuffer> cluster_prefixes(new InputBuffer(prefix + prefix_suf, 
				all_intervals.find(prefix_suf)->second, buffer_size, 0) );
			prefixes.push_back(cluster_prefixes);
			shared_ptr<InputBuffer> cluster_suffixes(new InputBuffer(prefix + suffix_suf,
				all_intervals.find(suffix_suf)->second, buffer_size, 0) );
			suffixes.push_back(cluster_suffixes);
			i++;
		}
		seen_so_far.resize(cores.size() + 1, 0);

		string other_suf = ".quals.other.lz";
		// cerr << "onto the general: " << (prefix + other_suf) << endl;
		// cerr << (all_intervals.find(other_suf) != all_intervals.end() ) << endl;
		other_qvs = shared_ptr<InputBuffer>(new InputBuffer(prefix + other_suf, 
			all_intervals.find(other_suf)->second, buffer_size, 0) );
	}

	// overloaded base function
	pair<int, unsigned long> seekToBlockStart(int const ref_id, 
		int const start_coord, int const end_coord) {
		// seek to block start on the membership stream
		bool t = false;
		auto start = data_in->loadOverlappingBlock(ref_id, start_coord, end_coord, t);
		if (start.first < 0) {
			cerr << "[ERROR] Could not navigate to the begining of the interval" << endl;
			exit(1);
		}
		// now seek for all of our core, prefix, suffix streams
		for (auto &str : cores) str->loadOverlappingBlock(ref_id,start_coord,end_coord, t);
		for (auto &str : prefixes) str->loadOverlappingBlock(ref_id,start_coord,end_coord, t);
		for (auto &str : suffixes) str->loadOverlappingBlock(ref_id,start_coord,end_coord, t);

		other_qvs->loadOverlappingBlock(ref_id, start_coord, end_coord, t);
		return start;
	}

	////////////////////////////////////////////////////////
	string getNextQualVector() {
		// read from the membership vector
		int i = getNextMember();
		string q_v = "*";
		if (i == GENERIC_PILE_ID) {// the general pile
			seen_so_far[0]++;
			q_v = getLine(other_qvs);
		}
		else {
			seen_so_far[i]++;
			auto core = getLine(cores[i-1]);
			auto pref = getLine(prefixes[i-1]);
			auto suf = 	getLine(suffixes[i-1]);
			q_v = pref + explodeString(core) + suf;
		}
		return q_v;
	}
};

#endif