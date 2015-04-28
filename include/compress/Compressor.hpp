#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
// bit vector compression
// #include <sdsl/bit_vectors.hpp>

// #include <lzip.h>
// #include <compress.h>

#include "IOLibParser.hpp"
#include "IOLibAlignment.hpp"
#include "OutputBuffer.hpp"
#include "QualityCompressor.hpp"
#include "RefereeUtils.hpp"




////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
struct Output_args {
	shared_ptr<OutputBuffer> offsets_buf;
	shared_ptr<OutputBuffer> edits_buf;
	shared_ptr<OutputBuffer> left_clips_buf;
	shared_ptr<OutputBuffer> right_clips_buf;
	shared_ptr<OutputBuffer> ids_buf;
	shared_ptr<OutputBuffer> flags_buf;
	// shared_ptr<OutputBuffer> quals_buf;
	shared_ptr<QualityCompressor> quals_buf;
	shared_ptr<OutputBuffer> opt_buf;
	shared_ptr<OutputBuffer> unaligned_buf;

	void flush() {
		if (flags_buf != nullptr)
			flags_buf->flush();
		offsets_buf->flush();
		edits_buf->flush();
		left_clips_buf->flush();
		right_clips_buf->flush();
		if (ids_buf != nullptr)
			ids_buf->flush();
		if (quals_buf != nullptr)
			quals_buf->flush(); // causes individual clusters to flush their OutputBuffers; notify the courier
		if (opt_buf != nullptr)
			opt_buf->flush();
		unaligned_buf->flush();
	}
};



class Compressor {
public:
	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	Compressor (string const & file_name, int t, Output_args & output_buffers, bool aligned_seq_only, bool unique_seq_only):
		parser(file_name, t),
		out_buffers(output_buffers),
		file_name(file_name),
		aligned_seq_only(aligned_seq_only),
		unique_seq_only(unique_seq_only) {
		failed_ = parser.failed();
		// TODO: check that successfully created a file parser
	}

	bool failed() {return failed_;}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	void compress() {
		bool first = true;

		ofstream head_out(file_name + ".head");
		/*SAM_hdr*/ auto* h = parser.header();
		auto num_ref = h->nref;
		for (auto i = 0; i < num_ref; i++) {
			head_out << i << " " << h->ref[i].name << endl;
		}
		head_out.close();

		int line_id = 0;
	    while ( parser.read_next() ) {
	        bam_seq_t* read = parser.getRead();
	        IOLibAlignment al(read);
	        if ( al.isUnalined() ) {
	        	unaligned_cnt++;
	            processUnalignedRead(al);
	        }
	        else {
	            processRead(al, first);
	            first = false;
	        }
	        line_id++;
	    }
	    flushUnalignedReads();
	    parser.close();
	    cerr << "Of them unaligned: " << unaligned_cnt << endl;
	    if (unique_seq_only)
	    	cerr << "Unique total reads: " << count << endl;

	    outputPair(offset_pair);
	    processHasEditsBits(edit_flags, file_name);
	}

////////////////////////////////////////////////////////////////
//
// Private methods and variables
//
////////////////////////////////////////////////////////////////
private:

	bool failed_ = false;

	bool aligned_seq_only = false; // compress all aligned sequence, including the multimaps

	bool unique_seq_only = false;	// omit multimaps, record data for a given read only once
	// TODO: strategies for choosing the alignement: smallest errors, most consistent offsets

	int count = 0;

	// parser
	IOLibParser parser;
	int unaligned_cnt = 0;
	int prev_offset = 0;
	int prev_ref = -1;
	// records offset from the previous alignment and the number of times we observe this offset
	pair<int,int> offset_pair;

	vector<bool> edit_flags;

	// vector<IOLibAlignment> unaligned_reads;
	vector<UnalignedRead> unaligned_reads;

	Output_args out_buffers;

	string file_name;

	// store the mapping for the flags, mapq and rnext
	// TODO: write out at the end
	unordered_map<short,int> flags_map;
	unordered_map<int,int> mapq_map;
	unordered_map<int,int> rnext_map;
	int flags_i = 0, mapq_i = 0, rnext_i = 0;

	////////////////////////////////////////////////////////////////
	void writeLeftClip(IOLibAlignment & al) {
		auto clip_length = al.lsc();
		auto clip_bytes = al.getLeftSoftClip();
		// for 2 bit encoding
		// lc_stream << clip_length; 
		// for (auto i = 0; i < clip_length; i++) lc_stream << clip_bytes[i];
		// lc_stream << endl;
		writeClip(clip_bytes, clip_length, out_buffers.left_clips_buf);
	}

	////////////////////////////////////////////////////////////////
	void writeRightClip(IOLibAlignment & al) {
		auto clip_length = al.rsc();
		auto clip_bytes = al.getRightSoftClip();
		// for 2 bit encoding
		// rc_stream << clip_length;
		// for (auto i = 0; i < clip_length; i++) rc_stream << clip_bytes[i];
		// rc_stream << endl;
		writeClip(clip_bytes, clip_length, out_buffers.right_clips_buf);
	}

	////////////////////////////////////////////////////////////////
	// write out edits to a compressed stream
	////////////////////////////////////////////////////////////////
	bool handleEdits(IOLibAlignment & al) {
		if ( !al.isPrimary() && unique_seq_only ) return false;

		bool hasEdits = al.handleEdits();

		// cerr << "after hadnled edits: rej=" << al.isRejected() << " hasEd=" << hasEdits << endl;
		// return early if no edits or can not encode the edits for some reason
		if (al.isRejected() || !hasEdits) {
			return hasEdits;
		}
		// if (!hasEdits) {
		// 	addUnsignedByte(0, out_buffers.edits_buf);
		// 	return hasEdits;
		// }

	    int num_edit_bytes = 0;
	    // calculate number of bytes needed for edits
	    if (al.lsc() > 0) {
	    	// cerr << "L";
	    	num_edit_bytes++;
	    }
	    if (al.rsc() > 0) {
	    	// cerr << "R";
	    	num_edit_bytes++;
	    }
	    if (al.lhc() > 0) {
	    	// cerr << "l";
	    	num_edit_bytes += 2;
	    }
	    if (al.rhc() > 0) {
	    	// cerr << "r";
	    	num_edit_bytes += 2;
	    }
	    if (al.merged_edits.size() > 0) {
	    	// cerr << "checking merged edits... ";
	        for (auto pr : al.merged_edits) {
	            if (pr.edit_op == 'E') { 
	            	// cerr << "E";
	            	if (pr.splice_len > 65535 )
	            		num_edit_bytes += 5;
	            	else
	                	num_edit_bytes += 4;
	            }
	            else if (pr.edit_op != 'R' && pr.edit_op != 'L') {
	            	// cerr << pr.edit_op;
	                num_edit_bytes += 2;
	            }
	        }
	    }
	    // assert(num_edit_bytes < 256);
		if (num_edit_bytes > 256) {
	    	al.set_rejected(true);
	    	return hasEdits;
		}
		if (num_edit_bytes == 0) {// TODO: this should be caught up stream
			// write "no edits", move on to the next alignment
			// addUnsignedByte(num_edit_bytes, out_buffers.edits_buf);
			cerr << "no edits after all" << endl;
			return false;
		}
		// cerr << " Edit len: " << (int)num_edit_bytes << endl;
	    // edit_stream << (unsigned char)num_edit_bytes;
	    addUnsignedByte(num_edit_bytes, out_buffers.edits_buf);

	    // write out clipping events
	    int prev_edit_offset = 0;
	    if (al.lsc() > 0) {
	        // edit_stream << 'L';
	        prev_edit_offset += al.lsc();
	        addUnsignedByte('L', out_buffers.edits_buf);
	        writeLeftClip(al);
	    }
	    if (al.rsc() > 0) {
	        // edit_stream << 'R';
	        addUnsignedByte('R', out_buffers.edits_buf);
	        writeRightClip(al);
	    }
	    // if (al.lhc() > 0) edit_stream << 'l' << (unsigned char)al.lhc();
	    // if (al.rhc() > 0) edit_stream << 'r' << (unsigned char)al.rhc();
	    if (al.lhc() > 0) {
		    addUnsignedByte('l', out_buffers.edits_buf);
		    addUnsignedByte(al.lhc(), out_buffers.edits_buf);
		    prev_edit_offset += al.lhc();
		}
		if (al.rhc() > 0) {
		    addUnsignedByte('r', out_buffers.edits_buf);
		    addUnsignedByte(al.rhc(), out_buffers.edits_buf);
		}

		// cerr << al.offset() << " ";
	    // write out all other edits
	    if (al.merged_edits.size() > 0) {
	        for (auto edit : al.merged_edits) {
	        	// cerr << (char) edit.edit_op << "(" << edit.edit_pos << ")";
	            if (edit.edit_op != 'E') {
	            	if (edit.edit_op != 'R' && edit.edit_op != 'L') {// these are handled separately
	            		writeOp(edit, prev_edit_offset, out_buffers.edits_buf);
	                	// edit_stream << (unsigned char) edit.edit_op; // edit code            
	                	// edit_stream << (unsigned char) (edit.edit_pos - prev_edit_offset);    // record how many bases since the last edit
	                	prev_edit_offset = edit.edit_pos;
	            	}
	            }
	            else {
	                // write out 'E', offset into the read
	                int splice_len = edit.splice_len;
	                // al.splices.erase(al.splices.begin());
	                // edit_stream << edit.edit_op << (unsigned char) (edit.edit_pos - prev_edit_offset);
	                // write out length of intron
	                // writeShort(edit_stream, splice_len);
	                // cerr << "Splice offset: " << (int)(edit.edit_pos - prev_edit_offset) << " len: " << splice_len;
	                if (splice_len > 65535 ) {
                		// cerr << " (long splice)" << endl;
                		writeLongSpliceOp(edit, prev_edit_offset, splice_len, out_buffers.edits_buf);
	                }
	                else {
	                	// cerr << " (short splice)" << endl;
	                	writeSpliceOp(edit, prev_edit_offset, splice_len, out_buffers.edits_buf);
	                }
	                prev_edit_offset = edit.edit_pos;
	            }
	        }
	    }
	    // cerr << endl;
		return hasEdits;
	}

	void outputPair(pair<int,int> const & offset_pair) {
		if (offset_pair.second > 1) {// offset with a multiplier
    		addOffsetPair(offset_pair.first, offset_pair.second, out_buffers.offsets_buf);
		}
		else { // offset w/o a multiplier
			addOffset(offset_pair.first, out_buffers.offsets_buf);
		}
	}

	////////////////////////////////////////////////////////////////
	// 1 2 3 4
	// 1 (1,1) 1
	// - (1,2) 2
	// - (1,3) 3
	// - (1,4) 4
	// bug 1: last pair is lost
	// bug 2: can accumulate starting with the first 1
	////////////////////////////////////////////////////////////////
	void handleOffsets(IOLibAlignment & al, bool hasEdits, bool first) {
		// if only want to encode a single read once ever, will skip all secondary alignements
		if ( !al.isPrimary() && unique_seq_only ) return;
		count++;
		int offset = al.offset();
	    // same offset as before -- increment the counter
	    if (first) {
	    	offset_pair = make_pair(offset, 1);
	    }
	    else {
	    	auto delta = (offset - prev_offset);
	    	if (offset_pair.first == delta ) {
		    	offset_pair.second++;
		    	// cerr << "branch 2" << endl;
		    }
		    else {
		    	// offset is different -- write prev to byte stream, store this new one
		    	outputPair(offset_pair);
		    	// reset the counter
		    	offset_pair = make_pair(delta, 1);
		    }
	    }
	    prev_offset = offset;
	    edit_flags.push_back(hasEdits);
	}

	////////////////////////////////////////////////////////////////
	// transform the flags before writing out
	////////////////////////////////////////////////////////////////
	void handleFlags(IOLibAlignment & al) {
		// remap flags, mapq, and rnext to a smaller domain
		short flags = al.flags();
		int mapq = al.mapq(), rnext = al.rnext();
		if ( flags_map.find(flags) == flags_map.end() ) {
			flags_map[flags] = flags_i++; 	// insert and increment
		}
		flags = flags_map[flags];
		// new mapq value
		if ( mapq_map.find(mapq) == mapq_map.end() ) {
			mapq_map[mapq] = mapq_i++; 	// insert and increment
		}
		mapq = mapq_map[mapq];
		// new rnext value
		if ( rnext_map.find(rnext) == rnext_map.end() ) {
			rnext_map[flags] = rnext_i++; 	// insert and increment
		}
		rnext = rnext_map[rnext];
		writeFlags(flags, mapq, rnext, al.pnext(), al.pnext() - al.tlen(), out_buffers.flags_buf);
	}

	////////////////////////////////////////////////////////////////
	// unordered_map<string, int> dictionary;
	// size_t unique_chunk_id = 0;
	// string delimiters = ".:_\s#";
	void handleReadNames(IOLibAlignment & al) {
		// TODO: need to transform the read IDs before passing it down to the lzip
		writeName(al.read_name(), out_buffers.ids_buf);

		// split, diff, minimize!
		// this can make compression for read ids faster and provide marginal improvements
		// string read_name = al.read_name();
		// size_t left = 0;
		// size_t right = read_name.find_first_of(delimiter);
		// vector<int> remapped_name;
		// while (right != string::npos ) {
		// 	string chunk = read_name.substr(left, right);
		// 	// if have not seen this before -- add to dictionary, increment unique_id
		// 	if (dictionary.find(chunk) == dictionary.end() ) {
		// 		dictionary[chunk] = unique_chunk_id;
		// 		unique_chunk_id++;
		// 	}
		// 	// output the mapped id instead of the chunk itself
		// 	remapped_name.push_back(dictionary[chunk]);
		// 	left = read_name.find_first_not_of(delimiter, right + 1);
		// 	right = read_name.find_first_of(delimiter, left + 1);
		// }
		// writeName2(remapped_name, out_buffers.ids_buf);
	}
	// TODO: write out the dictionary

	// TODO: can only initialize parameters w/in constructors?
	// QualityCompressor qual_compressor;
	////////////////////////////////////////////////////////////////
	void handleQuals(IOLibAlignment & al) {
		if (al.isPrimary()) {
			// haven't seen the read before -- output the quals
			auto len = al.read_len();
			auto* q = al.quals();
			// flip the qual vector based on the reverse-complement bit
			if (al.isRC()) {
				for (int i = 0; i < len / 2; i++) {
					char temp = q[i];
					q[i] = q[len - i - 1];
					q[len - i - 1] = temp;
				}
			}
			out_buffers.quals_buf->handleRead(q, len);
			// qual_compressor.handleRead(q, len);

			// just plzip
			// writeQualVector(q, len, out_buffers.quals_buf);
		}
	}

	////////////////////////////////////////////////////////////////
	void handleOptionalFields(IOLibAlignment & al) {
		// opt_stream << al.opt_fields() << endl;
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	bool printed_warning = false;
	void processUnalignedRead(IOLibAlignment & al) {
		// save for later
		unaligned_reads.emplace_back(al.read_name(), al.read_name_len(), al.getSeq(), al.quals());
		// TODO: 10K - arbitrary parameter, may be as big or as small as one wants
		if (unaligned_reads.size() >= 10000) {
			flushUnalignedReads();
		}
	}

	void flushUnalignedReads() {
		// cerr << "flushing unaligned n=" << unaligned_reads.size() << endl;
		//for (auto r : unaligned_reads) { for (auto c : r.seq) cerr << c; cerr << endl;}
		sort(unaligned_reads.begin(), unaligned_reads.end(), [] (UnalignedRead const & a, UnalignedRead const & b) {
			// returns ​true if the first argument is less than (i.e. is ordered before) the second.
	//		for (auto c : b.seq) cerr << (char)c; cerr << " --" << endl;
	//		assert(a.seq.size() == b.seq.size());
			for (int i = 0; (i < a.seq.size() ) && (i < b.seq.size() ); i++) {
				if (a.seq[i] < b.seq[i]) return true;
				else if (a.seq[i] > b.seq[i]) return false;
			}
			// reads identical -- return false as per strict weak ordering
			return false;
		});
		// write out
		for (auto read : unaligned_reads) {
			writeUnaligned(read, out_buffers.unaligned_buf);
		}
		unaligned_reads.clear();
	}

	void processHasEditsBits(vector<bool> & has_edits, string & prefix) {
		// cerr << "TODO: save the has_edits bit vector" << endl;
	    // writeAsRRR(has_edits, prefix);
	    writeAsBytes(has_edits, prefix);
	}

	void writeAsBytes(vector<bool> & has_edits, string & prefix) {
		size_t size = has_edits.size();
		ofstream has_edits_out(prefix + ".has_edits");
		/* write the number of bits we expect -- may not need this since we 
		 have # alignments implicitly from offsets file */
		uint8_t mask = 255;
		has_edits_out << (uint8_t)( (size >> 24) & mask);
		has_edits_out << (uint8_t)( (size >> 16) & mask);
		has_edits_out << (uint8_t)( (size >> 8) & mask);
		has_edits_out << (uint8_t)(size & mask);
		// write out binary data as bytes
		int num_bytes = (int)ceil(size / 8.0);
		for (auto i = 0; i < num_bytes; i++) {
			uint8_t byte = 0;
			for (auto j = 0; (j < 8) && (i*8+j < size); j++) {
				byte = (byte << 1) | has_edits[i*8 + j];
			}
			has_edits_out << byte;
		}
		has_edits_out.close();
	}

	/*
	void writeAsRRR(vector<bool> & has_edits, string & prefix) {
		uint64_t size = has_edits.size();
	    // create bv of size bits
	    sdsl::bit_vector bv(size);

	    // set sparse bits
	    for(uint64_t i = 0; i < size; i++) {
	        bv[i] = has_edits[i];
	    }
	    //
	    // sdsl::rrr_vector<> rrr_vec63(bv);
	    // sdsl::rrr_vector<> rrr_vec31(bv);
	    sdsl::rrr_vector<> rrr_vec15(bv);
	    // sdsl::rrr_vector<> rrr_vec127(bv);

	    // std::cout << "bv uses " << sdsl::size_in_mega_bytes(bv) << " MB." << std::endl;
	    std::cout << "rrr_vec15 uses " << sdsl::size_in_mega_bytes(rrr_vec15) << " MB." << std::endl;
	    sdsl::store_to_file(rrr_vec15, prefix + ".has_edits");
	}
	*/

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	bool processRead(IOLibAlignment & al, bool & first) {
	    bool rejected = false, new_transcript = false;
	    auto ref = al.ref();
	    if (first) {
	    	cerr << "Assuming uniform read len of " << al.read_len() << " bases" << endl;
	    	writeReadLen(al.read_len(), out_buffers.edits_buf);
	    	// edit_stream << (uint8_t)al.read_len();
	    	prev_ref = ref;

	    	// write the reference id before writing out alignment offsets
	    	addReference(prev_ref, first, out_buffers.offsets_buf);
	    	new_transcript = true;
	    }
	    else if (ref != prev_ref) {
	    	// starting a different chromosome -- dump offsets to a file
	    	outputPair(offset_pair);
	    	offset_pair = make_pair(0,0);
	    	prev_offset = 0;

	    	addReference(ref, first, out_buffers.offsets_buf);
	    	prev_ref = ref;
	    	new_transcript = true;
	    }

	    if (!aligned_seq_only && !unique_seq_only)
	    	handleReadNames(al);
	    bool hasEdits = handleEdits(al);
	    if (al.isRejected()) {
	    	return true;
	    }
	    handleOffsets(al, hasEdits, first || new_transcript);
	    if (!aligned_seq_only && !unique_seq_only) {
	    	handleFlags(al);
	    	handleQuals(al);
	    	handleOptionalFields(al);
		}

	    first = false;
	    
	    return al.isRejected();
	}
};

#endif
