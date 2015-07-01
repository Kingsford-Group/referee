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
#include "ReadIDStream.hpp"
#include "FlagsStream.hpp"
#include "QualityStream.hpp"
// #include "MergedEditsStream.hpp"
#include "TranscriptsStream.hpp"



////////////////////////////////////////////////////////////////

#define D_SEQ 				1
#define D_FLAGS				2
#define D_READIDS			4
#define	D_QUALS				8
#define D_OPTIONAL_FIELDS	16


////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
struct InputStreams {
	shared_ptr<ClipStream> left_clips;
	shared_ptr<ClipStream> right_clips;
	shared_ptr<OffsetsStream> offs;
	shared_ptr<EditsStream> edits;

	shared_ptr<FlagsStream> flags;
	shared_ptr<ReadIDStream> readIDs;
	shared_ptr<QualityStream> qualities;

	InputStreams() {}
};

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
class Decompressor {

	unordered_map<int,string> t_map;

	////////////////////////////////////////////////////////////////////////////
	// sync the streams
	void sync_streams(InputStreams & is, pair<int, unsigned long> & off_block_start,
		pair<int, unsigned long> & edit_block_start, int target_coord) {
		int offset_coord = off_block_start.first;
		unsigned long offset_num_al = off_block_start.second;

		int edit_coord = edit_block_start.first;
		unsigned long edit_num_al = edit_block_start.second;
		cerr << "of " << offset_num_al << " e " << edit_num_al << endl;

		// first sync be the number of coordinates preceeding the blocks
		cerr << "bringing edits up to speed" << endl;
		while (edit_num_al < offset_num_al) {
			is.edits->next(); // next has_edit byte
			// skip edit list for this alignment if there is one
			if (is.edits->hasEdits() ) is.edits->getEdits();
			edit_num_al++;
		}
		cerr << "edit_num_al: " << edit_num_al << endl;
		while (offset_num_al < edit_num_al) {
			// need off to seek forward to edit_start
			is.offs->getNextOffset();
			offset_num_al++;
		}
		cerr << "offs_num_al: " << offset_num_al << endl;

		// now seek to the target coordinate
		int offset = is.offs->getCurrentOffset();
		assert(offset_coord <= offset);
		while (offset < target_coord) {
			offset = is.offs->getNextOffset();
			if (offset == END_OF_TRANS) {
				cerr << "[ERROR] unexpected end of transcript offsets" << endl;
				exit(1);
			}
			else if (offset == END_OF_STREAM) {
				cerr << "[ERROR] unexpected end of offsets stream" << endl;
				exit(1);
			}
			// advance to the next edit
			is.edits->next();
			if (is.edits->hasEdits() ) is.edits->getEdits();
		}
		cerr << "after seeking current coord is: " << offset << endl;
	}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
public:

	Decompressor(
		string const & input_fname,
		string const & output_fname,
		string const & ref_path,
		unordered_map<int,string> const & Ts):
		file_name(input_fname),
		output_name(output_fname),
		ref_path(ref_path),
		t_map(Ts) { }

	////////////////////////////////////////////////////////////////////////////
	// Reconstruct SAM file by combining the inputs; restoring reads and quals
	////////////////////////////////////////////////////////////////////////////
	void decompress(int read_len, InputStreams & is, uint8_t const options) {
		// sequence-specific streams
		// read_len = is.edits->getReadLen();
		cerr << "[decompress the entire contents]" << endl;
		cerr << "Read length:\t" << (int)read_len << endl;
		assert(read_len > 0);
		TranscriptsStream transcripts(file_name, ref_path, "-d", t_map);
		recovered_file.open( output_name.c_str() );
		check_file_open(recovered_file, file_name);

		// prime the first blocks in every stream
		is.offs->seekToBlockStart(-1, 0, 0);
		is.edits->seekToBlockStart(-1, 0, 0);
		is.left_clips->seekToBlockStart(-1, 0, 0);
		is.right_clips->seekToBlockStart(-1, 0, 0);
		if (is.flags == nullptr) cerr << "Flags stream is not set up" << endl;
		is.flags->seekToBlockStart(-1, 0, 0);

		int ref_id = is.offs->getCurrentTranscript();
		int i = 0;
		while ( is.offs->hasMoreOffsets() ) {
			// zero-based offsets
			int offset_0 = is.offs->getNextOffset();
			// cerr << "Offset_1: " << (offset_0 + 1);

			if (offset_0 == END_OF_TRANS) {
				// remove the prev transcript sequence -- will not need it anymore
				transcripts.dropTranscriptSequence(ref_id);
				// pull the next transcript -- seq will be loaded on first access
				ref_id = is.offs->getNextTranscript();
				if (ref_id == END_OF_STREAM) {
					cerr << "Done" << endl;
					return;
				}
				cerr << "chr=" << transcripts.getMapping(ref_id) << " ";
			}
			else if (offset_0 == END_OF_STREAM) {
				// break
				cerr << "done";
			}
			else {
				// legit offset
				int ret = is.edits->next(); // advance to the next alignment
				if (ret == END_OF_STREAM) {
					cerr << "no more edits" << endl;
				}

				reconstructAlignment(offset_0, read_len, ref_id, transcripts,
					is.edits,
					is.left_clips, is.right_clips, is.readIDs, is.flags, is.qualities,
					options);
			}
			i++;
			if (i % 1000000 == 0) {
				cerr << i / 1000000 << "mln ";
			}
		}
		cerr << endl;
		recovered_file.close();
	}


	////////////////////////////////////////////////////////////////
	// Decompress alignments within a given interval
	////////////////////////////////////////////////////////////////
	void decompressInterval(GenomicInterval interval, int read_len, InputStreams & is,
		const uint8_t options) {
		TranscriptsStream transcripts(file_name, ref_path, "-d", t_map);
		recovered_file.open( output_name.c_str() );
		if (!recovered_file) {
			cerr << "[ERROR] Could not open output file." << endl;
			exit(1);
		}
		cerr << "Read len=" << read_len << endl;
		assert(read_len > 0);
		cerr << "[decompress alignments from an interval]" << endl;

		interval.chromosome = transcripts.getID(to_string(interval.chromosome) );
		int ref_id = interval.chromosome;
		// TODO: can return false when no data for that interval is available
		pair<int,unsigned long> off_start_coord = is.offs->seekToBlockStart(interval.chromosome, interval.start, interval.stop);
		assert(off_start_coord.first == is.offs->getCurrentOffset() );
		// loads the first block overlapping he requested coordinate
		// cerr << "Seeking to the first edit block" << endl;
		pair<int,unsigned long> edit_start_coord = is.edits->seekToBlockStart(interval.chromosome, interval.start, interval.stop);

		// cerr << "syncing input streams" << endl;
		is.left_clips->seekToBlockStart(-1, 0, 0);
		is.right_clips->seekToBlockStart(-1, 0, 0);

		sync_streams(is, off_start_coord, edit_start_coord, interval.start);
		// cerr << "STREAMS SYNCED" << endl;

		// now restore alignments
		int i = 0, offset = 0;
		while ( is.offs->hasMoreOffsets() ) {
			offset = is.offs->getNextOffset();
			// cerr << "offset: " << offset << endl;

			if (offset >= interval.stop) {
				cerr << "[INFO] Reached the end of the interval (" << interval.stop << " <= " << offset << ")" << endl;
				return;
			}

			if (offset == END_OF_TRANS) {
				ref_id = is.offs->getNextTranscript();
				if (ref_id == END_OF_STREAM) {
					cerr << "Done" << endl;
					return;
				}
				// cerr << "chr=" << transcripts.getMapping(ref_id) << " ";
			}
			else if (offset == END_OF_STREAM) {
				// break
				cerr << "done";
			}
			else {
				// legit offset
				int ret = is.edits->next(); // advance to the next alignment
				if (ret == END_OF_STREAM) {
					// cerr << "done with edits" << endl;
					// break;
				}
				// cerr << "has edit: " << edits.hasEdits() << endl;

				reconstructAlignment(offset, read_len, ref_id, transcripts,
					is.edits,
					is.left_clips, is.right_clips, is.readIDs, is.flags, is.qualities,
					options);

			}
			i++;
			if (i % 1000000 == 0) {
				cerr << i / 1000000 << "mln ";
			}
		}
		cerr << endl;
		recovered_file.close();
	}

////////////////////////////////////////////////////////////////
//
// Privates
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
	// offset is 0-based
	////////////////////////////////////////////////////////////////
	void reconstructAlignment(int offset, int read_len, int ref_id,
			TranscriptsStream & transcripts,
			shared_ptr<EditsStream> edits,
			shared_ptr<ClipStream> left_clips,
			shared_ptr<ClipStream> right_clips,
			shared_ptr<ReadIDStream> read_ids,
			shared_ptr<FlagsStream> flags,
			shared_ptr<QualityStream> qualities,
			uint8_t const options) {

		string cigar, md_string, read;
		bool has_edits = edits->hasEdits();
		if (has_edits) {
			md_string = "MD:Z:";
			vector<uint8_t> edit_ops = edits->getEdits();
			read = buildEditStrings(read_len, edit_ops, cigar, md_string,
				left_clips, right_clips, offset, ref_id, transcripts);
		}
		else {	
			read = transcripts.getTranscriptSequence(ref_id, offset, read_len);
		}
		// to upper case
		// std::transform(read.begin(), read.end(), read.begin(), ::toupper);


		if (options & D_READIDS) {
			string read_id;
			if ( read_ids->getNextID(read_id) != SUCCESS) {
				// cerr << "[INFO] Could not extract the next read ID. Using the default." << endl;
				read_id = "*";
			}
			recovered_file << read_id << "\t";
		}

		int flag = -1, mapq = -1, rnext = -1, pnext = -1, tlen = -1;
		if (options & D_FLAGS) {
			if (flags != nullptr) {
				auto alignment_flags = flags->getNextFlagSet();
				assert(alignment_flags.size() == 5);
				flag = alignment_flags[0]; 
				mapq = alignment_flags[1];
				rnext = alignment_flags[2]; 
				pnext = alignment_flags[3],
				tlen = alignment_flags[4];
			}

			recovered_file << flag << "\t";
			// write out reference name, offset (SAM files use 1-based offsets)
			recovered_file << transcripts.getMapping(ref_id) << "\t" << (offset + 1);
			recovered_file << "\t" << mapq;
			if (has_edits)
				recovered_file << "\t" << cigar;
			else
				recovered_file << "\t" << (int)read_len << "M";
			recovered_file << "\t";
			if (rnext < 0) {
				recovered_file << "*";
				pnext = 0;
				tlen = 0;
			}
			else if (rnext == 0) {
				recovered_file << "=";
				tlen = pnext - tlen;
			}
			else recovered_file << rnext;

			recovered_file << "\t" << pnext << "\t" << tlen << "\t";
		}

		if (options & D_SEQ) {
			recovered_file << read;
				// more data to come -- separate
			if ( (options & D_QUALS) || (options & D_OPTIONAL_FIELDS) )
				recovered_file << "\t";
		}
		// TODO: write out qual vector
		if (options & D_QUALS) {
			recovered_file << qualities->getNextQualVector();
		}
		if (options & D_OPTIONAL_FIELDS) {
			if (has_edits)
				cerr << md_string << " ";
			// TODO: write out other optional fields
		}
		recovered_file << endl;
	}

	////////////////////////////////////////////////////////////////
	// build CIGAR, MD strings for the read with edits
	// offset -- zero based
	////////////////////////////////////////////////////////////////
	string buildEditStrings(int read_len, vector<uint8_t> & edits,
		string & cigar, string & md_string,
		shared_ptr<ClipStream> left_clips,
		shared_ptr<ClipStream> right_clips,
		int offset, int ref_id,
		TranscriptsStream & transcripts) {

		string right_cigar, right_clip, read = transcripts.getTranscriptSequence(ref_id, offset, read_len);
		int j = 0, last_md_edit_pos = 0, last_cigar_edit_pos = 0, clipped_read_len = read_len;
		int offset_since_last_cigar = 0; // reset to 0 when on the cigar edit, incremented when on MD edit
		int offset_since_last_md = 0;
		int splice_offset = 0;
		int last_abs_pos = 0, Ds = 0, Is = 0; // number of deletions
		bool first_md_edit = true;
		
		// cerr << " len=" << edits.size() << endl;
		// if (offset == 57509325) {
		// 	for (auto e : edits)
		// 		cerr << (int)e << " ";
		// 	cerr << endl;
		// }

		uint8_t op;
		while (j < edits.size() ) {
			op = edits[j];
			// if (offset == 57509325)
				// cerr << op << " ";
			switch (op) {
				case 'L': {
					string left_clip;
					if (left_clips == nullptr) {
						// clipped data not available
						break;
					}
					else {
						left_clips->getNext(left_clip);
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
				}
				break;
				case 'R': {
					if (right_clips == nullptr) {
						// clipped data not available
						break;
					}
					else {
						right_clips->getNext(right_clip);
						right_cigar += to_string(right_clip.length());
						right_cigar += "S";
						clipped_read_len -= right_clip.length();
					}
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
					if (offset == 57509325) cerr << offset << endl;
					if (read.size() <= last_cigar_edit_pos) {
						cerr << "weird stuff: " << read.size() << " cigar: " << last_cigar_edit_pos << " offs: " << offset << endl;
					}
					assert(read.size() > last_cigar_edit_pos);
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
				assert(read_len - right_clip.length() < read.size() );
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

		// cerr << "Done w/ edits";
		return read;
	}

};

#endif
