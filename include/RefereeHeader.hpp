#ifndef REFEREE_HEADER_LIB
#define REFEREE_HEADER_LIB

class RefereeHeader {

	string version;

	// path to the head file
	string path;

	unordered_map<int,string> t_map;

	unordered_map<int,short> flags_map;

	unordered_map<int,int> mapq_map;

	unordered_map<int,int> rnext_map;

	// lengths in bases for each reference sequence (as recorded in *.head file)
	// TODO: fill out
	unordered_map<int,size_t> lengths;

	int read_len;

	pair<int,int> parseFlagLine(string const & line) {
		// cerr << line << endl;
		auto idx = line.find(" ");
		// cerr << idx << " ";
		auto idx2 = line.find(" ", idx+1);
		// cerr << idx2 << " ";
		int original = stoi(line.substr(idx+1, idx2 - idx) );
		// cerr << original << " ";
		int index = stoi(line.substr(idx2+1));
		// cerr << index;
		return make_pair(index, original);
	}


public:
	RefereeHeader(string const & fname) : path(fname) {}

	////////////////////////////////////////////////////////////////
	// TODO: parse using libstaden?
	// mapping from ref_id to a reference_name (e.g. 0 -> chr10)
	////////////////////////////////////////////////////////////////
	unordered_map<int, string> parse() {
		ifstream f_in(path);
		check_file_open(f_in, path);
		cerr << "[INFO] Reading reference sequence name mapping." << endl;

		int t_index = 0;
		string line, t_name, type, chromo;
		while (getline(f_in, line)) {
			if (line.find("HD") != string::npos) {
				// version
				auto idx = line.find(separator);
				version = line.substr(idx+1);
			}
			else if (line.find("read_len=") != string::npos ) {
				// read len
				auto idx = line.find("=");
				read_len = stoi(line.substr(idx+1));
			}
			else if (line.find("flags") != string::npos) {
				auto p = parseFlagLine(line);
				flags_map.insert(p);
			}
			else if (line.find("mapq") != string::npos) {
				auto p = parseFlagLine(line);
				mapq_map.insert(p);
			}
			else if (line.find("rnext") != string::npos) {
				auto p = parseFlagLine(line);
				rnext_map.insert(p);
			}
			else {
				auto idx = line.find_first_of(separator);
				type = line.substr(0, idx);
				t_index = stoi(type);
				auto idx2 = line.find_first_of(separator, idx+1);
				chromo = line.substr(idx+1, idx2-idx-1);
				t_map[t_index] = chromo;
				int length = stol(line.substr(idx2+1));
				lengths[t_index] = length;
			}
		}
		f_in.close();
		return t_map;
	}

	////////////////////////////////////////////////////////////////
	vector<int> getTranscriptIDs() {
		vector<int> all_ids;
		for (auto p : t_map) all_ids.push_back(p.first);
		sort(all_ids.begin(), all_ids.end());
		return all_ids;
	}

	////////////////////////////////////////////////////////////////
	string getMapping(int ref_id) {
		return t_map[ref_id];
	}

	string get_version() {return version;}

	int getReadLen() { return read_len;}

	size_t getTranscriptLength(int t_id) {
		if (lengths.find(t_id) == lengths.end()) return -1;
		return lengths[t_id];
	}

	unordered_map<int, string> getTranscriptIDsMap() {return t_map;}

	unordered_map<int,short> getFlagsEncoding() {return flags_map;}

	unordered_map<int,int> getMapqEncoding() { return mapq_map;}

	unordered_map<int,int> getRnextEncoding() { return rnext_map;}
};

#endif