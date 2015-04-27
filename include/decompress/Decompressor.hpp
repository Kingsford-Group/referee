#ifndef DECOMPRESSOR_HPP
#define DECOMPRESSOR_HPP

#include <unordered_map>
#include <queue>
// #include <memory>

#include <iostream>
#include <fstream>
#include <cassert>
#include <queue>

// #include "InputBuffer.hpp"
#include "OffsetsStream.hpp"
#include "EditsStream.hpp"
#include "ClipStream.hpp"
#include "MergedEditsStream.hpp"
#include "TranscriptsStream.hpp"


char reverseReplace(uint8_t & c) {
        if (c == 'V') return 'N';
        if (c == 'W') return 'A';
        if (c == 'X') return 'C';
        if (c == 'Y') return 'G';
        if (c == 'Z') return 'T';
        if (c == 'K') return 'A';
        if (c == 'L') return 'C';
        if (c == 'M') return 'G';
        if (c == 'N') return 'T';
        return '-';
}

class Decompressor {
////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
public:

	Decompressor(string const & input_fname, 
		string const & output_fname,
		string const & ref_path):
		file_name(input_fname),
		output_name(output_fname),
		ref_path(ref_path) { }

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	void decompress() {
		OffsetsStream offs(file_name);
		EditsStream edits(file_name);
		ClipStream left_clips(file_name, ".left_clip" );
		ClipStream right_clips(file_name, ".right_clip" );
		// MergedEditsStream edits(file_name, read_len);
		read_len = edits.getReadLen();
		cerr << "Read length:\t" << (int)read_len << endl;
		TranscriptsStream transcripts(file_name, ref_path);
		recovered_file.open( output_name.c_str() );
		// TODO: check if opened successfully

		int ref_id = offs.getNextTranscript();
		int i = 0;
		while ( offs.hasMoreOffsets() ) {
			int offset = offs.getNextOffset();
			if (offset == END_OF_TRANS) {
				ref_id = offs.getNextTranscript();
				if (ref_id == END_OF_STREAM) {
					cerr << "Done" << endl;
					return;
				}
				cerr << "chr=" << transcripts.getMapping(ref_id) << " ";
			}
			else if (offset == END_OF_STREAM) {
				// break
				// cerr << "done";
			}
			else {
				// legit offset
				int ret = edits.next(); // advance to the next alignment
				if (ret == END_OF_STREAM) {
					cerr << "done with edits" << endl;
					// break;
				}
				// cerr << "has edit: " << edits.hasEdits() << endl;
				if (edits.hasEdits() ) {
					// extract edits
					vector<uint8_t> edit_ops = edits.getEdits();
					reconstructRead(offset, read_len, ref_id, transcripts, edit_ops, left_clips, right_clips); // TODO: add right and left clips
				}
				else {
					reconstructRead(offset, ref_id, transcripts);
				}
			}
			i++;
			// if (i > 20) exit(1);
			if (i % 1000000 == 0) {
				cerr << i / 1000000 << "mln ";
			}
		}
		recovered_file.close();
	}

////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
private:

	string file_name; // file prefix for the files with parts of the data

	string output_name; // file prefix for the files with parts of the data

	string ref_path; // path to the file with the reference sequence

	uint8_t read_len; // uniform read length

	ofstream recovered_file;

	////////////////////////////////////////////////////////////////
	// reconstructs read without edits
	// TODO: fill out a iolib staden BAM record and write to a bam format
	////////////////////////////////////////////////////////////////
	void reconstructRead(int offset, int ref_id, TranscriptsStream & transcripts) {
		// SAM files use 1-based offsets
		recovered_file << transcripts.getMapping(ref_id) << "\t" << (offset+1);
		string read = transcripts.getTranscriptSequence(ref_id, offset, read_len);
		std::transform(read.begin(), read.end(), read.begin(), ::toupper);
		// CIGAR string is all matches
		recovered_file << "\t" << (int)read_len << "M";
		// recovered_file << "\t" << read; /*<< /*"\tMD:Z:" << (int)read_len;*/
		recovered_file << endl;
	}

	////////////////////////////////////////////////////////////////
	// reconstructs a read that had edits
	////////////////////////////////////////////////////////////////
	void reconstructRead(int offset, int read_len, int ref_id, TranscriptsStream & transcripts, vector<uint8_t> & edits, 
			ClipStream & left_clips,
			ClipStream & right_clips) {
		recovered_file << transcripts.getMapping(ref_id) << "\t" << (offset+1) ;
		// TODO
		// int left_clip_len = 0, right_clip_len = 0;
		// string read = getTranscriptSequence(ref_id, offset + left_clip_len, read_len - left_clip_len - right_clip_len);
		string cigar, md_string = "MD:Z:";
		string read = buildEditStrings(read_len, edits, cigar, md_string, left_clips, right_clips, offset, ref_id, transcripts);
		std::transform(read.begin(), read.end(), read.begin(), ::toupper);
		recovered_file << "\t" << cigar;
		// recovered_file << "\t" << read << /*"\t" << md_string <<*/ endl;
		recovered_file << endl;
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	string buildEditStrings(int read_len, vector<uint8_t> & edits, string & cigar, string & md_string, ClipStream & left_clips,
			ClipStream & right_clips,
			int offset, int ref_id, TranscriptsStream & transcripts) {
		string right_cigar, right_clip, read = transcripts.getTranscriptSequence(ref_id, offset, read_len);
		// cerr << "Read before edits: " << read << endl;
		int j = 0, last_md_edit_pos = 0, last_cigar_edit_pos = 0, clipped_read_len = read_len;
		int offset_since_last_cigar = 0; // reset to 0 when on the cigar edit, incremented when on MD edit
		int offset_since_last_md = 0;
		int splice_offset = 0;
		int last_abs_pos = 0, Ds = 0, Is = 0; // number of deletions
		bool first_md_edit = true;
		// cerr << endl << "off " << offset+1 << " len=" << edits.size() << " ";
		uint8_t op;
		while (j < edits.size() ) {
			op = edits[j];
			// cerr << op << " ";
			switch (op) {
				case 'L': {
					string left_clip;
					left_clips.getNext(left_clip);
					// just concatenate the left clip and the read
					read = left_clip + read.substr(0, read_len - left_clip.length());
					// update cigar string
					cigar += to_string(left_clip.length());
					cigar += "S";
					// update counters
					last_cigar_edit_pos += left_clip.length();
					last_abs_pos += left_clip.length();
					offset_since_last_cigar = 0;
					offset_since_last_md = 0;
					last_md_edit_pos = 0;
					clipped_read_len -= left_clip.length();
				}
				break;
				case 'R': {
					right_clips.getNext(right_clip);
					right_cigar += to_string(right_clip.length());
					right_cigar += "S";
					clipped_read_len -= right_clip.length();
				}
				break;
				case 'l': {
					j++;
					// update the read (shorten)
					read = read.substr(edits[j]);
					// update cigar string
					cigar += to_string(edits[j]);
					cigar += "H";
					last_cigar_edit_pos += edits[j];
					last_abs_pos += edits[j];
					offset_since_last_cigar = 0;
					offset_since_last_md = 0;
					clipped_read_len -= edits[j];
				}
				break;
				case 'r': {
					j++;
					cigar += to_string(edits[j]);
					cigar += "H";
					last_cigar_edit_pos += edits[j];
					clipped_read_len -= edits[j];
				}
				break;
				case 'D': {
					// handle a deletion
					int ds = 0;
					bool first_d = true;
					offset_since_last_cigar += edits[j+1] - Is;
					while (j < edits.size()) {// consume all Ds
						if (edits[j] != 'D') {
							// end of run of D's
							break;
						}
						else if (!first_d && edits[j+1] > 0 ) {
							// non contiguous deletions
							break;
						}
						ds++;
						// cerr << (char)edits[j] << "-" << (int)edits[j+1] << ",";
						last_abs_pos += edits[j+1];
						read.erase(last_abs_pos, 1);
						first_d = false;
						j += 2;
					}
					j--;
					assert(j < edits.size());
					Ds += ds;
					// cerr << "D's pos: " << last_abs_pos << " #=" << ds << " ";

					// compensate for deleted bases by adding to the end of the read from the reference
					// splice_offset: handles the case when this might be after a splicing event
					// cerr << "seq past the read: " << transcripts.getTranscriptSequence(ref_id, offset + splice_offset + read_len, 5) << " ";
					read += transcripts.getTranscriptSequence(ref_id, offset + splice_offset + read_len, ds);
					// cerr << read << " ";
					if (offset_since_last_cigar > 0) {
						cigar += to_string(offset_since_last_cigar);
						cigar += "M";
					}
					cigar += to_string(ds);
					cigar += "D";
					last_cigar_edit_pos += offset_since_last_cigar;
					offset_since_last_cigar = 0;
				}
					break;
				case 'E':
				case 197: {
					// handle a splice
					// cerr << "splice ";
					bool is_long_splice = edits[j] >> 7;
					j++;
					offset_since_last_cigar += edits[j] - Is;
					last_abs_pos += edits[j];
					// cerr << "at pos:" << last_abs_pos << " ";

					// update cigar string
					cigar += to_string(offset_since_last_cigar);
					cigar += "M";
					// get splice len
					// cerr << "len: ";
					j++;
					int splice_len = ( (int)edits[j] << 8);
					// cerr << "j=" << j << " " << (int)edits[j] << " ";
					j++;
					splice_len |= (int)edits[j];
					// cerr << "j=" << j << " " << (int)edits[j] << " ";
					if (is_long_splice) {
						j++;
						splice_len = splice_len << 8;
						splice_len |= (int)edits[j];
						// cerr << "j=" << j << " " << (int)edits[j] << " ";
					}
					cigar += to_string(splice_len);
					cigar += "N";
					// update the counters
					splice_offset += offset_since_last_cigar + splice_len;
					last_cigar_edit_pos += offset_since_last_cigar;
					offset_since_last_cigar = 0;
					// update the read
					// TODO: test with hard clipped strings
					read.replace(last_cigar_edit_pos, read_len - last_cigar_edit_pos,
						transcripts.getTranscriptSequence(ref_id, offset + splice_offset, read_len - last_cigar_edit_pos));
					// cerr << "read: " << read << endl;
				}
					break;
				case 'V':
				case 'W':
				case 'X':
				case 'Y':
				case 'Z': {
					/* insert bases into the reference to recover the original read */
					// TODO: if several Is in a row -- need to handle them together
					// if Is are disjoint -- will handle them separately
					// TODO: is more than one I in a row possible?
					int is = 0;
					bool first_i = true;
					offset_since_last_cigar += edits[j+1] - Is;
					while (j < edits.size()) {// consume all Is
						if (edits[j] < 'V' || edits[j] > 'Z') {
							// not an I -- break the loop
							// cerr << "break loop 1 ";
							break;
						}
						else if (!first_i && edits[j+1] > 0) {
							// non contiguous Is -- break the loop
							// cerr << "break loop 2 j=" << j << " ";
							break;
						}
						is++;
						// cerr << (char)edits[j] << "-" << (int)edits[j+1] << ",";
						last_abs_pos += edits[j+1];
						// insert the missing base
						read.insert(read.begin() + last_abs_pos, reverseReplace(edits[j]) );// insert a single char
						j += 2;
						first_i = false;
					}
					j--;
					assert(j < edits.size());
					// cerr << "I's pos: " << last_abs_pos << " ";
					// trim the read to be of appropriate length
					read.resize(read_len);
					// update the global insertion counter
					Is += is;

					if (offset_since_last_cigar > 0) {
						cigar += to_string(offset_since_last_cigar);
						cigar += "M";
					}
					cigar += to_string(is);
					cigar += "I";
					last_cigar_edit_pos += offset_since_last_cigar + is;
					offset_since_last_cigar = 0;
				}
				break;
				default: {
					// TODO: handle mismatches
					j++;
					if (first_md_edit) {
						// cerr << "(+" << (int)edits[j] + Is << ") ";
						md_string += to_string(edits[j]);
						first_md_edit = false;
						last_md_edit_pos += edits[j] + 1;
					}
					else {
						// cerr << "(+" << (int)edits[j]-1+Is << ") ";
						md_string += to_string(edits[j] - 1);
						last_md_edit_pos += edits[j];
					}
					last_abs_pos += edits[j] + Is;
					// add the letter we see in the reference to the MD string
					// cerr << "last abs pos: " << last_abs_pos << " ";
					md_string += read[last_abs_pos];
					// add correct read to have its original base 
					read[last_abs_pos] = op;
					offset_since_last_cigar += edits[j];
					offset_since_last_md = 0;
				}
			}
			j++;
		}
		// TODO: update read w/ right soft/hard clip if it took place
		if (last_cigar_edit_pos < read_len) {
			if (right_cigar.size() > 0) {
				cigar += to_string(read_len - last_cigar_edit_pos - right_clip.length());
				cigar += "M";
				cigar += right_cigar;
				// update read w/ the right clip
				read.replace(read_len - right_clip.length(), right_clip.length(), right_clip);
			}
			else {
				cigar += to_string(read_len - last_cigar_edit_pos);
				cigar += "M";
			}
		}
		if (last_md_edit_pos < clipped_read_len) {
			md_string += to_string(clipped_read_len - last_md_edit_pos);
		}
		// cerr << cigar << "\t" << md_string << endl;
		return read;
	}


	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	// string applyEdits(int const & offset, char const & read_len, 
	// 	string const & edit,
	// 	bool rc,
	// 	shared_ptr<queue<string>> & left_clips,
	// 	shared_ptr<queue<string>> & right_clips,
	// 	shared_ptr<string> & trans_seq, string & read_id) {
	// 	string read = trans_seq->substr(offset, read_len);

	// 	auto edit_len = edit.size();
	// 	if (edit_len == 0)
	// 		return read;

	// 	bool print = true;
	// 	// if (read_id.compare("SOLEXA5_133_6_91_4706_2537_0") == 0) print = true;

	// 	vector<pair<char,unsigned short>> edit_pairs;
	// 	unsigned char offset_into_read = 0;
	// 	int edit_pos = 0;
	// 	int del_adjust = 0;
	// 	int ins_adjust = 0;
	// 	int j = 0;
	// 	string left_clip;
	// 	if (edit[j] == 'L') {// has left clip
	// 		left_clip = left_clips->front();
	// 		if (print) cerr << "L " << left_clip << " ";
	// 		left_clips->pop();
	// 		read = left_clip + read.substr(0, read_len - left_clip.size());
	// 		assert(read.size() == read_len);
	// 		edit_pos += left_clip.size();
	// 		// cerr << "ed pos: " << edit_pos << " "; 
	// 		j++;
	// 	}
	// 	if (j < edit_len && edit[j] == 'R') {
	// 		string right_clip = right_clips->front();
	// 		if (print) cerr << "R " << right_clip << " ";
	// 		right_clips->pop();
	// 		read.replace(read_len - right_clip.size(), right_clip.size(), right_clip);
	// 		assert(read.size() == read_len);
	// 		j++;
	// 	}

	// 	// deal with splices
	// 	while (j < edit_len && edit[j] == 'E') {
	// 		// TODO
	// 		j += 4;
	// 		cerr << "splice done ";
	// 	}

	// 	while (j < edit_len) {
	// 		char edit_letter = edit[j++];
	// 		edit_pos += edit[j];
	// 		// cerr << "ed pos: " << edit_pos << " "; 
	// 		if (print) cerr << edit_letter << "," << (int)edit[j] << " ";

	// 		if (edit_letter == 'D') {
	// 			assert(trans_seq->size() > offset - 1 + read.size() + 1 + del_adjust);
	// 			del_adjust++;
	// 		}
	// 		else if (edit_letter == 'H') {
	// 			cerr << "Hard clip" << endl;
				
	// 		}
	// 		else if (edit_letter >= 'V' && edit_letter <= 'Z') {
	// 			read.insert(edit_pos, 1, reverseReplace(edit_letter));
	// 			// once we're done with all the edits, we will erase the right most chars that slid off
	// 		}
	// 		else if (edit_letter == 'E') {
	// 			cerr << "this should not happen ";
	// 			if (print) { cerr << read_len - offset_into_read << " ";
	// 			cerr << "t: " << trans_seq->substr(offset - 1 + edit[j], read_len) << " ";
	// 			}
	// 			string spliced_piece = trans_seq->substr(offset - 1 + offset_into_read - left_clip.size() + edit[j], read_len - offset_into_read);
	// 			if (print) cerr << spliced_piece;
	// 			read.replace(offset_into_read + ins_adjust + del_adjust, spliced_piece.size(), spliced_piece);
	// 			if (print) cerr << " r: " << read << endl;
	// 		}
	// 		else {
	// 			// cerr << (edit_pos + ins_adjust) << " (" << read.size() << ") ";
	// 			read[edit_pos + ins_adjust] = edit_letter;
	// 		}
	// 		j++;
	// 	}
	// 	if (print)cerr << endl;

	// 	// trim read at the end if we have accumulated additional bases through insertions
	// 	read = read.substr(0, read_len);
	// 	if (rc) {
	// 		reverse_complement(read);
	// 	}
	// 	return read;
	// }

};

#endif