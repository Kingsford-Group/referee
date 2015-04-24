#ifndef TRANS_STREAM_H
#define TRANS_STREAM_H

#include <unordered_map>

#include "TransUtils.hpp"

const string separator = "\t\s ";

class TranscriptsStream {
public:
	////////////////////////////////////////////////////////////////
	TranscriptsStream (string const & file_name, string const & ref): ref_path(ref) {
		t_map = parseTranscriptIDsPlain(file_name + ".head");
		readFAI(ref + ".fai");
		ref_sequence = parseTranscripts(ref_path);
	}

	////////////////////////////////////////////////////////////////
	string getMapping(int ref_id) {
		return t_map[ref_id];
	}

	////////////////////////////////////////////////////////////////
	// offset -- 0-based offset into the reference sequence
	// len -- lenght of the sequence to extract
	////////////////////////////////////////////////////////////////
	string getTranscriptSequence(int const ref_id, int const offset, int const len) {
		// TODO: use fai to read one seq at a time
		auto mapped_name = t_map[ref_id];
		auto it = ref_sequence.find(mapped_name);
		if (it == ref_sequence.end() ) {
			cerr << "[ERROR] Could not find sequence for " << mapped_name << endl;
			return "";
		}
		if (it->second->size() < offset + len) {
			cerr << "[ERROR] Offset is past the length of the reference sequence" << endl;
			return "";
		}
		return it->second->substr(offset, len);
	}

private:
	unordered_map<int, string> t_map;

	unordered_map<string,size_t> ref_offsets;

	// seuqence for the reference
	unordered_map<string, shared_ptr<string>> ref_sequence;

	// path to the reference
	string ref_path;

	////////////////////////////////////////////////////////////////
	// read index for the reference and stores offsets
	////////////////////////////////////////////////////////////////	
	void readFAI(string const & fname) {
		// TODO: fill out ref_offsets
	}

	////////////////////////////////////////////////////////////////
	// TODO: parse using libstaden
	// mapping from ref_id to a reference_name (e.g. 0 -> chr10)
	////////////////////////////////////////////////////////////////
	unordered_map<int, string> parseTranscriptIDs(string const & fname) {
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
			// cerr << line << endl;
			auto idx = line.find_first_of(separator);
			type = line.substr(0, idx);
			if (type.compare("@SQ") == 0) { // reference sequence line
				idx = line.find_first_not_of(separator, idx);
				auto idx2 = line.find_first_of(separator, idx);
				chromo = line.substr(idx, idx2 - idx);
				auto mid = chromo.find(":");
				t_name = chromo.substr(mid+1, idx2 - idx - mid);
				// t_map[t_name] = t_index;
				t_map[t_index] = t_name;
				cerr << "Reference " << t_name << "\t->\t" << t_index << endl;
				t_index++;
			}
		}

		f_in.close();
		return t_map;
	}

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
			// cerr << line << endl;
			auto idx = line.find_first_of(separator);
			type = line.substr(0, idx);
			t_index = stoi(type);
			chromo = line.substr(idx+1);
			t_map[t_index] = chromo;
			// cerr << "Reference " << chromo << "\t->\t" << t_index << endl;
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
};

#endif