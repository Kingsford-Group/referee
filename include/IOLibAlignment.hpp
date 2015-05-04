// Alignment class

#ifndef IO_LIB_ALIGNMENT_HPP
#define IO_LIB_ALIGNMENT_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <stdlib.h>

#include "RefereeUtils.hpp"
#include "TranscriptsStream.hpp"

using namespace std;

extern "C" {
    #include "io_lib/scram.h"
    #include "io_lib/os.h"
    #undef max
    #undef min
}


struct edit_pair {
	unsigned char edit_op = '\0';
	int splice_len = 0;
	int edit_pos = -1;

	edit_pair(unsigned char const & edit_op, int const & edit_pos) {
		this->edit_op = edit_op;
		this->edit_pos = edit_pos;
	}

	edit_pair(unsigned char const & edit_op, int const & edit_pos, int sl):
		splice_len(sl) {
		this->edit_op = edit_op;
		this->edit_pos = edit_pos;
	}
};

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
char mapIndels(char c) {
	if (c=='N') return 'V';
	if (c=='A') return 'W';
	if (c=='C') return 'X';
	if (c=='G') return 'Y';
	if (c=='T') return 'Z';
	return '-';
}

char twoBitToChar[] = {'A','C','G','T'};
char bit2char(uint8_t bits) {
	return twoBitToChar[samToTwoBit[bits]];
}

////////////////////////////////////////////////////////////////
// return true if the edit operation is a mismatch (one of A, C, G, T, N)
////////////////////////////////////////////////////////////////
inline bool isMismatch(unsigned char const & edit_op) {
    return edit_op == 'A' || edit_op == 'C' || edit_op == 'G' || edit_op == 'T' || edit_op == 'N';
}

struct UnalignedRead {
	vector<uint8_t> seq;
	vector<uint8_t> qual;
	vector<uint8_t> read_name;

	UnalignedRead(char * rn, int name_len, vector<uint8_t> s, char * q): seq(s) {
		for (auto i = 0; i < name_len; i++)
			read_name.push_back(rn[i]);
		auto q_len = s.size();
		for (auto i = 0; i < q_len; i++)
			qual.push_back(q[i] + '!');
	}
};

class IOLibAlignment {

////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
public:

	vector<edit_pair> merged_edits;

	// vector<int> splices;

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	IOLibAlignment(bam_seq_t* read): read(read) {

	}

	////////////////////////////////////////////////////////////////
	bool isPrimary() {
		return (bam_flag(read) & 0x900) == 0;
	}

	////////////////////////////////////////////////////////////////
	bool isRC() {
		return (bam_flag(read) & 16) == 0;
	}

	////////////////////////////////////////////////////////////////
	bool isUnalined() {
		return (bam_flag(read) & BAM_FUNMAP) > 0;
	}

	////////////////////////////////////////////////////////////////
	uint32_t ref() {return bam_ref(read); }

	////////////////////////////////////////////////////////////////
	uint32_t offset() { return bam_pos(read); }

	int read_name_len () const {return bam_name_len(read); }

	////////////////////////////////////////////////////////////////
	char* read_name () const { return bam_name(read); }

	////////////////////////////////////////////////////////////////
	int read_len() { return bam_seq_len(read); }

	////////////////////////////////////////////////////////////////
	short flags() { return bam_flag(read); }

	////////////////////////////////////////////////////////////////
	uint8_t mapq() { return bam_map_qual(read); }

	////////////////////////////////////////////////////////////////
	int rnext() { return bam_mate_ref(read); }

	////////////////////////////////////////////////////////////////
	int pnext() { return bam_mate_pos(read); }

	////////////////////////////////////////////////////////////////
	int tlen() { return bam_ins_size(read); }

	////////////////////////////////////////////////////////////////
	char * quals() { return bam_qual(read); }

	////////////////////////////////////////////////////////////////
	// w/o MD string
	char * opt_fields() { 
		auto auxillary_fields = bam_aux(read);
		// TODO
		// return bam_qual(read); 
		// return opt_fields;
		return auxillary_fields;
	}

	int lsc() { return left_soft_clip;}

	int rsc() { return right_soft_clip;}

	int lhc() { return left_hard_clip;}

	int rhc() { return right_hard_clip;}

	////////////////////////////////////////////////////////////////
	vector<uint8_t> getLeftSoftClip() {
		auto seq = bam_seq(read);
		vector<uint8_t> char_seq(left_soft_clip, 0);
		for (int i = 0; i < left_soft_clip; i++) {
			char_seq[i] = bit2char(bam_seqi(seq, i));
		}
		return char_seq;
	}

	////////////////////////////////////////////////////////////////
	// 2-bit encoding: <num_bytes> + char seq of 2-bit-encoded packed seq
	////////////////////////////////////////////////////////////////
	vector<uint8_t> getLeftSoftClip2() {
		auto seq = bam_seq(read);
		unsigned char num_bytes = (left_soft_clip >> 2) + ( (left_soft_clip & 3) > 0 );
		vector<uint8_t> bit_seq(num_bytes, 0);
		for (int i = 0; i < left_soft_clip; i++) {
			bit_seq[ i >> 2 ] = bit_seq[ i >> 2 ] | ( samToTwoBit[bam_seqi(seq, i)] << (2 * (i & 3) ) );
		}
		return bit_seq;
	}

	////////////////////////////////////////////////////////////////
	// sequence for right clip
	////////////////////////////////////////////////////////////////
	vector<unsigned char> getRightSoftClip() {
		auto seq = bam_seq(read);
		auto seq_len = bam_seq_len(read);
		vector<unsigned char> char_seq(right_soft_clip, 0);
		for (int i = 0; i < right_soft_clip; i++) {
			char_seq[i] = bit2char(bam_seqi(seq, seq_len - right_soft_clip + i));
		}
		return char_seq;
	}

	////////////////////////////////////////////////////////////////
	// 2-bit encoding: <num_bytes> + char seq of 2-bit-encoded packed seq
	////////////////////////////////////////////////////////////////
	vector<unsigned char> getRightSoftClip2() {
		auto seq = bam_seq(read);
		auto seq_len = bam_seq_len(read);
		auto num_bytes = (right_soft_clip >> 2) + ( (right_soft_clip & 3) > 0 );
		// char * bit_seq = malloc( num_bytes * sizeof(unsigned char) );
		vector<unsigned char> bit_seq(num_bytes, 0);
		for (int i = 0; i < right_soft_clip; i++) {
			bit_seq[ i >> 2 ] = bit_seq[ i >> 2 ] | ( samToTwoBit[bam_seqi(seq, seq_len - right_soft_clip + i)] << (2 * (i & 3) ) );
		}
		return bit_seq;
	}

	////////////////////////////////////////////////////////////////
	// get read's sequence
	////////////////////////////////////////////////////////////////
	vector<uint8_t> getSeq() const {
		auto seq = bam_seq(read);
		auto seq_len = bam_seq_len(read);
		vector<uint8_t> char_seq(seq_len, 0);
		for (int i = 0; i < seq_len; i++) {
			char_seq[i] = bit2char(bam_seqi(seq, i));
		}
		return char_seq;
	}

	bool isRejected() {return rejected;}

	void set_rejected(bool t) {rejected = t;}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	bool handleEdits(TranscriptsStream & ts) {
		auto cigarEdits = handleCigar(); // return smth?
	    auto mdEdits = handleMDstring(ts);
	    if (rejected) return cigarEdits || mdEdits;
	    // cerr << auxillary_fields << endl;
	    if (cigarEdits || mdEdits )
	    	merge();
	    return cigarEdits || mdEdits;
	}

private:
	bam_seq_t * read;

	unsigned char left_soft_clip = 0;
	unsigned char left_hard_clip = 0;

	unsigned char right_soft_clip = 0;
	unsigned char right_hard_clip = 0;

	bool rejected = false;

	// maybe replace w/ something else to save on mem allocs
	vector<edit_pair> md_edits;
	vector<edit_pair> cigar_edits;

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	bool handleCigar() { // return smth?
		auto cigar_len = bam_cigar_len(read);
	    auto cigar = bam_cigar(read);
	    // uint8_t* qseq = reinterpret_cast<uint8_t*>(bam_seq(read));
		if (cigar_len > 0) {
	    	size_t readIdx = 0;
	    	auto transcriptIdx = bam_pos(read);
	    	// if the read starts before the beginning of the transcript,
		    // only consider the part overlapping the transcript
		    if (transcriptIdx < 0) {
		        readIdx = -transcriptIdx;
		        transcriptIdx = 0;
		    }
	    	size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);

	    	uint32_t opOffset = 0;
	    	for (uint32_t cigarIdx = 0; cigarIdx < cigar_len; ++cigarIdx) {
	            uint32_t opLen = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
	            auto op = cigar[cigarIdx] & BAM_CIGAR_MASK;
	            consumeCigarOp(op, opLen, opOffset);
	        }
	    }
	    // if (cigar_edits.size() > 0) cerr << " (" << cigar_edits.size() << ")" << endl;
	    return cigar_edits.size() > 0;
	}



	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	void consumeCigarOp(int32_t op, uint32_t opLen, uint32_t & offset_into_read) {
	    switch (op) {
	        case BAM_UNKNOWN:
	            cerr << "[ERROR] Unknown symbol in the cigar string." << endl;
	            break;
	        case BAM_CMATCH:
	        case BAM_CBASE_MATCH:
	           // do nothing
	        	offset_into_read += opLen;
	            break;
	        case BAM_CBASE_MISMATCH:
	            // TODO
	        	offset_into_read += opLen;
	            break;
	        case BAM_CINS: {
	        		// cerr << "I ";
		        	int ins_offset = 0;
		        	for (auto i = 0; i < opLen; i++) {
						// assert(this->seq.size() > read_pos);
						auto seq = bam_seq(read);
						auto base_insert = bit2char(bam_seqi(seq, offset_into_read) );
						cigar_edits.push_back( edit_pair(mapIndels(base_insert), offset_into_read - ins_offset) );
						offset_into_read++;
						ins_offset++;
					}
				}
	            break;
	        case BAM_CDEL:
				for (auto i = 0; i < opLen; i++) {
					// cerr << "D ";
					cigar_edits.push_back( edit_pair('D', offset_into_read) );
				}
	            break;
	        case BAM_CREF_SKIP: {
	        	// the splicing event
	        	// cerr << "E ";
	        	cigar_edits.push_back( edit_pair('E', offset_into_read, opLen) );
	        	// cerr << "splice: " << opLen << "; ";
	        	// splices.push_back(opLen);
	        }
	        break;
	        case BAM_CSOFT_CLIP:
	        	if (offset_into_read > 0) { // S is the last character in the cigar string
	        		// cerr << "S ";
					right_soft_clip = opLen; // this->seq.substr(read_pos);
					if (right_soft_clip > 256) rejected = true;
					cigar_edits.push_back( edit_pair('R', opLen) );
				}
				else {
					// cerr << "S ";
					offset_into_read+= opLen;
					left_soft_clip = opLen;
					if (left_soft_clip > 256) rejected = true;
					cigar_edits.push_back( edit_pair('L', opLen) );
				}
	        	break;

	        case BAM_CHARD_CLIP:
	        	if (offset_into_read > 0) {// not on the right side
	        		// cerr << "H ";
					right_hard_clip = opLen;
					cigar_edits.push_back( edit_pair('r', opLen) );
				}
				else {
					// cerr << "H ";
					left_hard_clip = opLen;
					cigar_edits.push_back( edit_pair('l', opLen) );
				}
	            offset_into_read += opLen;
	            break;
	        case BAM_CPAD: // ???
	            break;
	    }
	}

	
    ////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	bool handleMDstring(TranscriptsStream & ts) {
		bool hasEdits = false;
		auto auxillary_fields = bam_aux(read);
		auto md = bam_aux_find(read, "MD");
		if (md != NULL) {
			md++; // skip 'Z' -- indicator of the printable string
			parseMD(md);
		}
		else {
			// check for mismatches
			md_edits = getMismatches(ts, cigar_edits);
		}
		hasEdits = md_edits.size() > 0;
		return hasEdits;
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	vector<edit_pair> getMismatches(TranscriptsStream & ts, vector<edit_pair> & cigar_edits) {
		vector<edit_pair> md_edits;
		vector<uint8_t> read_seq = this->getSeq();
		string ref_seq = ts.getTranscriptSequence(this->ref(), this->offset(), read_seq.size());
		for (int i = 0; i < read_seq.size(); i++) {
			if (read_seq[i] != ref_seq[i]) {
				// mismatch -- record the letter that appear in the read at position i
				md_edits.push_back( edit_pair(read_seq[i], i) ); 
			}
		}
		// TODO: handle clipped regions
		// TODO: handle splicing events in the cigar
		// TODO: handle indels in the cigar
		return md_edits;
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	void parseMD(const char * md) {
		auto mdLen = strlen(md);
		int read_pos = left_soft_clip;
		int offset = 0;
		for (auto i = 0; i < mdLen; i++) {
	        if ( isdigit(md[i]) ) {
	        	offset = offset * 10 + (md[i] - '0');
	        }
	        else {
	        	read_pos += offset;
	        	consumeMDOp(md, i, read_pos);
				offset = 0;
	        }
	    }
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	void consumeMDOp(const char * md, int & i, int & read_pos) {
		auto* seq = bam_seq(read);
		switch (md[i]) {
			case '^':
	    		// this the same as D in CIGAR: is this a duplication of edits?
				// "deletion from the reference" -- need to delete letters from the ref to get the read
				// e.g. gap in the read relative to the reference
				i++;
				while ( isalpha(md[i]) ) i++;
			break;
			case 'A': case 'C': case 'G': case 'T': case 'N': {
				// mismatch between the reference and the read: get the character in the read
				// that mismatched with the reference
				auto ch = bit2char(bam_seqi(seq, read_pos) );
				md_edits.push_back( edit_pair(ch, read_pos) ); // 0-based position for this mismatch in the clipped read
				read_pos++;
				// cerr << "(" << md_edits.back().edit_pos << "," << md_edits.back().edit_pos << ") ";
			}
			break;
			default:
				cerr << "[ERROR] Unexpected symbol (" << md[i] << ") in the MD string" << endl;
				md_edits.clear();
				return;
    	}
	}

	

    ////////////////////////////////////////////////////////////////
    //
    ////////////////////////////////////////////////////////////////
    void merge() {
    	if (cigar_edits.size() > 0 || md_edits.size() > 0) {
    		// spot_checked++;
	  //   	cerr << read_name() << " ";
	  //   	cerr << "cigar: ";
	  //   	for (auto p : cigar_edits) cerr << "(" << p.edit_op << "," << p.edit_pos << ") ";
			// cerr << "md: ";
			// for (auto p : md_edits) cerr << "(" << p.edit_op << "," << p.edit_pos << ") ";
	  //   	cerr << endl;
			if (cigar_edits.size() > 0) {
				if (md_edits.size() > 0) {
					std::merge(	cigar_edits.begin(), cigar_edits.end(), 
							md_edits.begin(), md_edits.end(), 
							std::back_inserter(merged_edits),
							[](edit_pair const & a, edit_pair const & b) {
								return a.edit_pos < b.edit_pos;
							});
					// cigar_edits.clear();
					// md_edits.clear();
				}
				else {
					merged_edits = cigar_edits;
					// md_edits.clear();
				}
			}
			else {
				merged_edits = md_edits;
				// cigar_edits.clear();
			}

			// some adjustments for insertion/deletion
			int ins = 0;
			// int clipped_bases = this->left_clip.size() + this->left_hard_clip;
			// cerr << "merged: ";
			auto seq = bam_seq(read);
			for (auto p : merged_edits) {
				if (p.edit_op >= 'V' && p.edit_op <= 'Z') ins++;
				if (isMismatch(p.edit_op) && ins > 0) {
					// cerr << "ins adj ";
					p.edit_op = bit2char( bam_seqi(seq, p.edit_pos + ins) );
					// p.edit_pos += ins;
					// cerr << "(" << p.edit_op << "," << p.edit_pos << ") ";
				}
			}
			// cerr << endl;
		}
	}
};

#endif




