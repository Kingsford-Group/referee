#ifndef TRANS_STREAM_H
#define TRANS_STREAM_H

#include <unordered_map>
#include <sstream>
#include <list>

#include "RefereeUtils.hpp"
#include "FastaReader.h"

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

	unordered_map<string, int> reverse_map;

	unordered_map<string,FaiEntry> fai_index;

	unordered_map<string,size_t> ref_offsets;

	// seuqence for the reference
	unordered_map<string, shared_ptr<string>> ref_sequence;

	// path to the reference
	string ref_path;

	int read_len = 0;

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	void indexFile(string const & fname, unordered_map<string,FaiEntry> & map) {
		ifstream f_in(fname);
		check_file_open(f_in, fname);
		// get file length
		f_in.seekg(0, f_in.end);
		size_t flen = f_in.tellg();
		f_in.seekg(0, f_in.beg);

		int bp_per_line = 0;
		int bytes_per_line = 0;
		ofstream fai_out(fname + ".fai");

		int block_size = pow(2,20);
		vector<uint8_t> buf( block_size );
		size_t bytes_read = 0;

		int64_t num_bases = 0, byte_offset = 0;
		string ref_name = "";
		// read all bytes in the file
		while (bytes_read < flen) {
			f_in.read( (char *) &buf[0], block_size);
			auto actually_read = f_in.gcount();
			bytes_read += actually_read;
			int i = 0;
			while (i < actually_read) {
				// new transcript
				if (buf[i] == '>') {
					if (i > 1) {
						// this is not the first first >
						map.insert( make_pair(ref_name,
							FaiEntry(ref_name, num_bases, byte_offset, bp_per_line, bytes_per_line) ) );
						// write out details of the previous trancsript
						fai_out << ref_name << "\t" << num_bases << "\t" << byte_offset << "\t" <<
							bp_per_line << "\t" << bytes_per_line << endl;
					}
					// find end of line, get ref_name
					i++; // skip '>'
					ref_name = "";
					while (buf[i] != '\n') {
						if (buf[i] == ' ') break;
						ref_name.push_back(buf[i]);
						i++;
					}
					i++; // skip end of line
					// past the header line -- here is where base sequence starts
					byte_offset = f_in.tellg();
					byte_offset -= (block_size - i);
					num_bases = 0;
					bp_per_line = 0;
					bytes_per_line = 0;
					while (buf[i] != '\n') {
						if ( isalpha(buf[i]) ) {
							bp_per_line++;
							num_bases++;
						}
						i++;
						bytes_per_line++;
					}
					bytes_per_line++;
					// done estimating bases and bytes per line
 				}
				else {
					if ( isalpha(buf[i]) ) {
						num_bases++;
					}
				}
				i++;
			}
			buf.clear();
		}
		// write out data for the last transcript
		int lines = (int) ceil ( (double) num_bases / bp_per_line );
		int newline_chars = bytes_per_line - bp_per_line;
		byte_offset = flen - (num_bases + newline_chars * lines);
		map.insert( make_pair(ref_name,
				FaiEntry(ref_name, num_bases, byte_offset, bp_per_line, bytes_per_line)
			) );
		fai_out << ref_name << "\t" << num_bases << "\t" << byte_offset << "\t" <<
			bp_per_line << "\t" << bytes_per_line << endl;

		fai_out.close();
	}

	////////////////////////////////////////////////////////////////
	// read index for the reference and stores offsets
	////////////////////////////////////////////////////////////////
	unordered_map<string,FaiEntry> readFAI(string const & fname) {
		unordered_map<string,FaiEntry> fai_index;
		ifstream fai_in(fname + ".fai");
		if (!fai_in) {
			cerr << "[INFO] " << fname << " .fai reference index not found. Creating one..." << endl;
			fai_index.clear();
			indexFile(fname, fai_index);
			cerr << "[INFO] Index written to " << fname << ".fai" << endl;
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
				fai_index.insert( make_pair( ref_name,
						FaiEntry(ref_name, n, b, bp, bytes)
						) );
			}
			fai_in.close();
		}
		cerr << "Found " << fai_index.size() << " index entries" << endl;
		return fai_index;
	}

	////////////////////////////////////////////////////////////////
	// TODO: parse using libstaden
	// mapping from ref_id to a reference_name (e.g. 0 -> chr10)
	////////////////////////////////////////////////////////////////
	unordered_map<int, string> parseTranscriptIDsPlain(string const & fname) {
		unordered_map<int,string> t_map;
		reverse_map.clear();
		ifstream f_in(fname);
		check_file_open(f_in, fname);
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
				reverse_map[chromo] = t_index;
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
	shared_ptr<string> readTranscriptSequence(string const & ref_name, FaiEntry & entry) {
		ifstream f_in(ref_name);
		if (!f_in) {
			cerr << "[ERROR] Could not open reference sequence file: " << ref_name << endl;
			exit(1);
		}
		int num_lines = (int) ceil( (double)entry.num_bases / entry.bases_per_line );
		int newline_chars = entry.bytes_per_line - entry.bases_per_line;
		int bytes_to_read = entry.num_bases + num_lines * newline_chars;
		// cerr << ref_name << ": bases " << entry.num_bases << " bytes to read " << bytes_to_read << endl;
		assert( entry.num_bases <= bytes_to_read);
		vector<char> S(bytes_to_read, 0);			// initialize a vector that is long enough
		f_in.seekg(entry.byte_offset);				// seek to the first base
		f_in.read( (char*)&S[0], bytes_to_read);	// read all bytes representing the sequence

		// remove newline characters
		auto start_time = chrono::system_clock::now();
		int i = 0, prev_i = 0, N = bytes_to_read;
		string str;
		while (i < N) {
			if ( isalpha(S[i]) ) {
				i++;
			}
			else {
				str.append(S.begin() + prev_i, S.begin() + i);
				i += newline_chars;
				prev_i = i;
			}
		}
		auto end_time = chrono::system_clock::now();
		assert(str.size() == entry.num_bases);
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
		fai_index = readFAI(ref);
	}

	////////////////////////////////////////////////////////////////
	string getMapping(int ref_id) {
		return t_map[ref_id];
	}

	int getID(string const & mapped_name) {
		return reverse_map[mapped_name];
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
			cerr << "[INFO] Loading sequence for " << mapped_name;
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
			cerr << " - loaded." << endl;
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
