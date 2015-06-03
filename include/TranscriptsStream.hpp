#ifndef TRANS_STREAM_H
#define TRANS_STREAM_H

#include <unordered_map>
#include <sstream>
#include <list>

const string separator = "\t\s ";

class TranscriptsStream {

	////////////////////////////////////////////////////////////////
	class FaiEntry {
	public:
		string ref_name;
		int64_t num_bases;
		int64_t byte_offset;
		int bases_per_line;
		int bytes_per_line;

		FaiEntry() {}

		FaiEntry(string s, int64_t n, int64_t b, int bp, int bytes):
			ref_name(s),
			num_bases(n),
			byte_offset(b),
			bases_per_line(bp),
			bytes_per_line(bytes) {}
	};

	unordered_map<int, string> t_map;

	unordered_map<string,FaiEntry> fai_index;

	unordered_map<string,size_t> ref_offsets;

	// seuqence for the reference
	unordered_map<string, shared_ptr<string>> ref_sequence;

	// path to the reference
	string ref_path;

	int read_len = 0;

	////////////////////////////////////////////////////////////////
	// read index for the reference and stores offsets
	////////////////////////////////////////////////////////////////	
	unordered_map<string,FaiEntry> readFAI(string const & fname) {
		unordered_map<string,FaiEntry> fai_index;
		ifstream fai_in(fname);
		if (!fai_in) {
			cerr << "[INFO] " << fname << " reference index not found. Creating one..." << endl;
			// TODO: fill out ref_offsets
		}
		else {
			// read it in
			cerr << "[INFO] Reading reference sequence index..." << endl;
			// this is usually a small file -- let's use stringstream for simplicity
			string line;
			
			while (getline(fai_in, line)) {
				istringstream ss(line);
				// have 5 tokens
				string ref_name, token;
				getline(ss, ref_name, '\t');
				// length in bases of this transcript
				getline(ss, token, '\t');
				auto n = stoul(token);
				// byte offset from the begining of the file to the first base of the sequence
				getline(ss, token, '\t');
				auto b = stoul(token);
				// bases per line
				getline(ss, token, '\t');
				int bp = stoi(token);
				// bytes per line
				getline(ss, token, '\t');
				int bytes = stoi(token);
				// fai_index.emplace(
				// 	std::piecewise_construct,
				// 	// std::forward_as_tuple(ref_name), 
				// 	// std::forward_as_tuple(ref_name, n, b, bp, bytes) );
				// 	ref_name,
				// 	FaiEntry{ref_name, n, b, bp, bytes} );
				fai_index.insert(
					make_pair(
						ref_name,
						FaiEntry(ref_name, n, b, bp, bytes)
						) );
			}
			cerr << "Fai index size: " << fai_index.size() << endl;
		}
		return fai_index;
	}

	////////////////////////////////////////////////////////////////
	// TODO: parse using libstaden
	// mapping from ref_id to a reference_name (e.g. 0 -> chr10)
	////////////////////////////////////////////////////////////////
	// unordered_map<int, string> parseTranscriptIDs(string const & fname) {
	// 	unordered_map<int,string> t_map;
	// 	ifstream f_in(fname);
	// 	if (!f_in) {
	// 		cerr << "[ERROR] Reference file not found: " << fname << endl;
	// 		exit(1);
	// 	}
	// 	cerr << "[INFO] Reading reference sequence name mapping." << endl; 

	// 	int t_index = 0;
	// 	string line, t_name, type, chromo;
	// 	while (getline(f_in, line)) {
	// 		// cerr << line << endl;
	// 		auto idx = line.find_first_of(separator);
	// 		type = line.substr(0, idx);
	// 		if (type.compare("@SQ") == 0) { // reference sequence line
	// 			idx = line.find_first_not_of(separator, idx);
	// 			auto idx2 = line.find_first_of(separator, idx);
	// 			chromo = line.substr(idx, idx2 - idx);
	// 			auto mid = chromo.find(":");
	// 			t_name = chromo.substr(mid+1, idx2 - idx - mid);
	// 			// t_map[t_name] = t_index;
	// 			t_map[t_index] = t_name;
	// 			cerr << "Reference " << t_name << "\t->\t" << t_index << endl;
	// 			t_index++;
	// 		}
	// 	}

	// 	f_in.close();
	// 	return t_map;
	// }

	////////////////////////////////////////////////////////////////
	// TODO: parse using libstaden
	// mapping from ref_id to a reference_name (e.g. 0 -> chr10)
	////////////////////////////////////////////////////////////////
	unordered_map<int, string> parseTranscriptIDsPlain(string const & fname) {
		unordered_map<int,string> t_map;
		ifstream f_in(fname);
		if (!f_in) {
			cerr << "[ERROR] Reference file not found: " << fname << endl;
			exit(1);
		}
		cerr << "[INFO] Reading reference sequence name mapping." << endl; 

		int t_index = 0;
		string line, t_name, type, chromo;
		while (getline(f_in, line)) {
			if (line.find("read_len=") == string::npos ) {
				auto idx = line.find_first_of(separator);
				type = line.substr(0, idx);
				t_index = stoi(type);
				chromo = line.substr(idx+1);
				t_map[t_index] = chromo;
				// cerr << "Reference " << chromo << "\t->\t" << t_index << endl;
			}
			else {
				// read len
				auto idx = line.find("=");
				read_len = stoi(line.substr(idx+1));
			}
		}

		f_in.close();
		return t_map;
	}


	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	inline unordered_map<string, shared_ptr<string>> parseTranscripts(string const & fname) {
		// check if file exists
		ifstream f_in(fname);
		if (!f_in) {
			cerr << "[ERROR] Reference file not found: " << fname << endl;
			exit(1);
		}
		f_in.close();

		unordered_map<string, shared_ptr<string>> transcripts;
		char * fname_char = new char[fname.size() + 1];
		strcpy(fname_char, fname.c_str());
		auto fr = FastaReader(fname_char);
		kseq_t * seq;
		while ((seq = fr.nextSequence())) {
			cerr << "[INFO] Reading reference sequence " << seq->name.s << "..." << endl;
			transcripts.emplace(seq->name.s, make_shared<string>(seq->seq.s) );
		}
		delete fname_char;

		return transcripts;
	}

	////////////////////////////////////////////////////////////////
	shared_ptr<string> readTranscriptSequence(string const & ref_name, FaiEntry & entry) {
		ifstream f_in(ref_name);
		if (!f_in) {
			cerr << "[ERROR] Could not open reference sequence file: " << ref_name << endl;
			exit(1);
		}
		int num_lines = (int) ceil( (double)entry.num_bases / entry.bases_per_line );
		int newline_chars = entry.bytes_per_line - entry.bases_per_line;
		int bytes_to_read = entry.num_bases + num_lines * newline_chars;
		cerr << ref_name << ": bases " << entry.num_bases << " bytes to read " << bytes_to_read << endl;
		assert( entry.num_bases <= bytes_to_read);
		vector<char> S(bytes_to_read, 0);			// initialize a vector that is long enough
		f_in.seekg(entry.byte_offset);				// seek to the first base
		f_in.read( (char*)&S[0], bytes_to_read);	// read all bytes representing the sequence

		// TODO: remove newline characters
		auto start_time = chrono::system_clock::now();
		// int i = 0, prev_i = 0, N = bytes_to_read;
		// string str;
		// while (i < N) {
		// 	if ( isalpha(S[i]) ) {
		// 		i++;
		// 	}
		// 	else {
		// 		str.append(S.begin() + prev_i, S.begin() + i);
		// 		i += newline_chars;
		// 		prev_i = i;
		// 	}
		// }
		int i = bytes_to_read - 1, prev_i = i;
		string str;
		while (i >= 0) {
			if ( isalpha(S[i]) ) {
				i--;
			}
			else {
				for (auto it = S.begin() + prev_i; it != S.begin() + i; i--)
					str.push_back(c);
				// pop these chars off
				i -= newline_chars;
				prev_i = i;
			}
		}
		reverse(str.begin(), str.end() );
		auto end_time = chrono::system_clock::now();
		cerr << "Trimming newlines: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << " us" << endl;
		cerr << "Sequence len: " << str.size() << endl;
		// append to the prev line
		return make_shared<string>(str);
	}

public:

	////////////////////////////////////////////////////////////////
	TranscriptsStream (string const & file_name, string const & ref, string const mode): ref_path(ref) {
		if (mode.compare("-c") == 0) {
			// if compressing: build a t_map
		}
		else {
			// if decompressing: read a t_map
			t_map = parseTranscriptIDsPlain(file_name + ".head");	
			// TODO: also: read dictionaries for the flags
		}
		fai_index = readFAI(ref + ".fai");
		// reads all sequences into memory at once
		// ref_sequence = parseTranscripts(ref_path);
	}

	////////////////////////////////////////////////////////////////
	string getMapping(int ref_id) {
		return t_map[ref_id];
	}

	////////////////////////////////////////////////////////////////
	void setMapping(int ref_id, string ref_name) {
		t_map[ref_id] = ref_name;
	}

	////////////////////////////////////////////////////////////////
	int getReadLength() {
		return read_len;
	}

	////////////////////////////////////////////////////////////////
	void dropTranscriptSequence(int const ref_id) {
		auto mapped_name = t_map[ref_id];
		if (ref_sequence.find(mapped_name) != ref_sequence.end() ) {
			ref_sequence.erase(mapped_name);
		}
	}

	////////////////////////////////////////////////////////////////
	// offset -- 0-based offset into the reference sequence
	// len -- lenght of the sequence to extract
	////////////////////////////////////////////////////////////////
	string getTranscriptSequence(int const ref_id, int const offset, int const len) {
		// use fai to read one seq at a time
		auto mapped_name = t_map[ref_id];
		auto it = ref_sequence.find(mapped_name);
		if (it == ref_sequence.end() ) {
			cerr << "[INFO] Could not find sequence for " << mapped_name << ". Loading... ";
			if ( fai_index.find(mapped_name) == fai_index.end() ) {
				cerr << "[ERROR] Reference name " << mapped_name << " not in the index." << endl;
				exit(1);
			}
			auto entry = fai_index[mapped_name];
			auto seq = readTranscriptSequence(ref_path, entry);
			if (seq->size() < offset + len) {
				cerr << "seq size: " << seq->size() << " vs  offs " << offset << " len " << len << endl;
			}
			assert(offset + len <= seq->size() );
			cerr << "Loaded." << endl;
			// store for fast access later
			ref_sequence[mapped_name] = seq;
			return seq->substr(offset, len);
		}
		if (it->second->size() < offset + len) {
			cerr << "[ERROR] Offset is past the length of the reference sequence" << endl;
			return "";
		}
		return it->second->substr(offset, len);
	}
};

#endif