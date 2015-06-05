#ifndef INPUT_STREAM_LIB_H
#define INPUT_STREAM_LIB_H

#include <vector>
#include <queue>
#include <deque>
#include <memory>
#include <fstream>
#include <system_error>
#include <map>

#include <fcntl.h>
#include <unistd.h>

#include <lzlib.h>

// #include "plzip/file_index.h"

// #include "tbb/concurrent_queue.h"

#include "IntervalTree.h"

using namespace std;
// using namespace tbb;

#define END_OF_STREAM -1
#define SUCCESS 1
#define END_OF_TRANS -2


// chromosome/transcript ID -- an interger
typedef int chromo_id_t;

#define chromo_min 0
#define chromo_max 300000000


////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
class MyBlock {
public:
	int compressed_size = -1;			// block's size (including header and trailer)
	int decompressed_size = -1;	// expected size of the decompressed data
	int offset = -1;			// offset to the block's header from the begining of the file

	MyBlock(int s, int exp, int off):
		compressed_size(s),
		decompressed_size(exp),
		offset(off) {}
};


////////////////////////////////////////////////////////////////
//
// http://www.nongnu.org/lzip/manual/lzlib_manual.html#Examples
// example 8
//
////////////////////////////////////////////////////////////////
vector<uint8_t> unzipData(shared_ptr<vector<uint8_t>> raw_data, int const new_data_size) {
	vector<uint8_t> unzipped_data(new_data_size, 0); // allocate needed amount of bytes

	// prepare a decoder stream
	LZ_Decoder * const decoder = LZ_decompress_open();
	if (LZ_decompress_errno( decoder ) != LZ_ok ) {
		LZ_decompress_close( decoder );
		cerr << "[ERROR] Could not initialize decoder" << endl;
		exit(1);
	}
	int wrote = 0, bytes_read = 0, need_to_write = raw_data->size();
	int chunk_size = LZ_decompress_write_size( decoder );
	assert(chunk_size > 0);
	chunk_size = std::min( chunk_size, need_to_write);

	bool trailing_garbage_found = false;
	while (wrote < need_to_write && !trailing_garbage_found && (chunk_size > 0) ) {
		int written = LZ_decompress_write( decoder, (uint8_t *) &(*raw_data)[wrote], chunk_size );
		if (written != chunk_size) {
			cerr << "[ERROR] Lzlib error: attempted to write (" << written << ")" << endl;
			exit(1);
		}
		if (written < 0) {
			cerr << "[ERROR] Could not write to decoder" << endl;
			exit(1);
		}
		wrote += written;
		int read = LZ_decompress_read( decoder, (uint8_t *) &unzipped_data[bytes_read], new_data_size - bytes_read );
		bytes_read += read;
		chunk_size = LZ_decompress_write_size( decoder );
	}
	assert( wrote >= need_to_write );
	assert(bytes_read == new_data_size);

	LZ_decompress_finish( decoder );
	auto retval = LZ_decompress_finished( decoder );
	// destroys all internal data
	LZ_decompress_close( decoder );
	return unzipped_data;
}


////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
class InputBuffer {

	bool no_blocks = false;

	// file stream name
	string name;

	map<chromo_id_t, IntervalTree<int,int> > chromosome_trees;

	int buffer_id;

	ifstream f_in;

	// a queue of blocks to go through and decompress one at a time
	deque<RawDataInterval> block_queue;

	// currently available decompressed bytes of the underlying data stream
	// (mostly consists of bytes from the most recently decompressed block)
	deque<uint8_t> bytes;

	int buffer_size;

	void readMoreLZIPBlocks() {
		if (block_queue.size() > 0) {
			auto block = block_queue.front();
			block_queue.pop_front();
			auto unzipped_data = decompressBlock(block);
			for (auto b : unzipped_data)
				bytes.push_back(b);
		}
		else {
			cerr << name << ": no more blocks" << endl;
		}
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	vector<MyBlock> seek_blocks(ifstream & f_in) {
		// cerr << "reading blocks from plziped stream" << endl;
		vector<MyBlock> blocks;
		// set cursor at the end of the stream
		f_in.seekg(0, f_in.end);
		// length of file
		int64_t file_len, pos, total = 0;
		file_len = pos = f_in.tellg();

		File_header header;
  		File_trailer trailer;
  		auto member_size = 0;
  		// scan until cursor reaches the begining of the file
  		while (pos > 0) {
  			f_in.seekg( -member_size - File_trailer::size, f_in.cur);
  			// read trailer
			f_in.read( (char *) trailer.data, File_trailer::size);
			// block size
			member_size = trailer.member_size();
			auto data_size = trailer.data_size();
			pos -= member_size;
			total += member_size;
			blocks.emplace_back(member_size, data_size, pos);
			try {
				// cerr << f_in.good() << " " << f_in.fail() << " " <<
					// f_in.bad() << " " << f_in.eof() << endl;
				f_in.exceptions(f_in.failbit);
			}
			catch (const std::ios_base::failure & e) {
				cerr << "[ERROR] Input error: " << e.what() << " " << /*e.code() <<*/ endl;
			}
		}
		assert(total == file_len);
		reverse( blocks.begin(), blocks.end() );
		return blocks;
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	vector<uint8_t> decompressBlock(RawDataInterval & block) {
		// read block bytes from the LZ stream
		// different non-C++11 compliant compilers may set failbit upon reaching eof
		// we reset all flags just in case
		f_in.clear();
		// set absolute position relative to the begining of the file
		f_in.seekg(block.byte_offset);
		try {
			f_in.exceptions(f_in.failbit);
		}
		catch (const std::ios_base::failure & e) {
			cerr << "[ERROR] Input error: " << e.what() << " " << /*e.code() <<*/ endl;
			cerr << "Try compiling against C++11-compliant standard libraries.";
			exit(1);
		}
		shared_ptr<vector<uint8_t>> raw_bytes(new vector<uint8_t>(block.block_size) );
		f_in.read( (char *) &(*raw_bytes)[0], block.block_size);
		auto unzipped_data = unzipData(raw_bytes, block.decompressed_size);
		cerr << name << ": decompressed " << block.decompressed_size << " bytes" << endl;
		return unzipped_data;
	}

	/////////////////////////////////////////////////////tellg2///////////
	//
	////////////////////////////////////////////////////////////////
	void createTree(int const chromo,
		vector<RawDataInterval> & chromo_intervals,
		map<chromo_id_t, IntervalTree<int,int> > & chromosome_trees) {
		chromosome_trees[chromo] =
				IntervalTree<int,int>(chromo_intervals);
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	void createChromosomeIntervalTree(
		shared_ptr<vector<TrueGenomicInterval>> genomic_intervals,
		vector<MyBlock> & lzip_blocks,
		map<chromo_id_t, IntervalTree<int,int> > & chromosome_trees) {

		assert(lzip_blocks.size() == genomic_intervals->size());

		int prev_chromo = genomic_intervals->at(0).start.chromosome;
		vector<RawDataInterval> chromo_intervals;
		for (int i = 0; i < genomic_intervals->size(); i++) {
			auto interval = genomic_intervals->at(i);
			auto block = lzip_blocks[i];
			if (prev_chromo != interval.start.chromosome) {
				// create a tree for intervals in [range_start, range_end]
				createTree(prev_chromo, chromo_intervals, chromosome_trees);
				// start a new range at this index
				chromo_intervals.clear();
				prev_chromo = interval.start.chromosome;
			}

			// interval.start.chromosome is prev_chromo
			if (interval.start.chromosome == interval.end.chromosome) {
				// start and end chromosome are the same and equal prev_chromo
				chromo_intervals.emplace_back(
					block.offset, block.compressed_size, block.decompressed_size,
					interval.start.chromosome, interval.start.offset, interval.end.offset, interval.num_alignments, interval.is_aligned);
				// moving on...
			}
			else {
				// figure out how many chromos this interval spans
				auto chromo_span = interval.end.chromosome - interval.start.chromosome - 1;
				// split the interval as many times as needed (at least 2 pieces)
				// first create a tree for chromo_intervals w/ a starting piece of the current interval
				chromo_intervals.emplace_back(
					block.offset, block.compressed_size, block.decompressed_size,
					interval.start.chromosome, interval.start.offset, chromo_max, interval.num_alignments, interval.is_aligned);
				// then create tree(s) for every chromo that is in between the start and the end
				createTree(interval.start.chromosome, chromo_intervals, chromosome_trees);
				chromo_intervals.clear();
				int k = 0;
				while (k < chromo_span) {
					chromo_intervals.emplace_back(
						block.offset, block.compressed_size, block.decompressed_size,
						interval.start.chromosome + k + 1, chromo_min, chromo_max, interval.num_alignments, interval.is_aligned);
					createTree(interval.start.chromosome + k + 1, chromo_intervals, chromosome_trees);
					chromo_intervals.clear();
					k++;
				}
				// then add the leftover piece to the chromo_intervals
				chromo_intervals.emplace_back(
					block.offset, block.compressed_size, block.decompressed_size,
					interval.end.chromosome, chromo_min, interval.end.offset, interval.num_alignments, interval.is_aligned);
				prev_chromo = interval.end.chromosome;
			}

		}
		// create a tree for intervals in chromo_intervals
		createTree(prev_chromo, chromo_intervals, chromosome_trees);
	}


////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
public:

	////////////////////////////////////////////////////////////////
	// default constructor: acquires file descriptor for a data stream,
	// builds a mapping between genomic coordinates and the byte blocks
	// within lzip-ed stream
	////////////////////////////////////////////////////////////////
	InputBuffer(string const & fname, shared_ptr<vector<TrueGenomicInterval>> genomic_intervals, int id, int bs = 1<<22):
		buffer_id(id),
		name(fname),
		buffer_size(bs),
		f_in(fname.c_str(), ifstream::in | ios::binary | ios::ate)  {
		check_file_open_silent(f_in, fname);
		// interval trees -- one per chromosome
		// fill out chromosome_trees
		auto lzip_blocks = seek_blocks(f_in);
		createChromosomeIntervalTree(genomic_intervals, lzip_blocks, chromosome_trees);
	}

	////////////////////////////////////////////////////////////////
	~InputBuffer() {
		f_in.close();
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	pair<int,unsigned long> loadOverlappingBlock(int const chromo, int const start_coord, int const end_coord,
		bool & is_transcript_start) {
		// cerr << "Loading an overlapping block for: " << name << endl;
		bytes.clear();
		block_queue.clear();

		if (chromo == -1) {
			// add all blocks from all tree (in order) to the block queue
			for (auto tree_p : chromosome_trees) {
				auto tree = tree_p.second;
				for (auto b : tree.intervals ) block_queue.push_back(b);
			}
			// decompress the first one
			RawDataInterval block = block_queue.front();
			is_transcript_start = block.isAlignedWithTranscriptStart();
			block_queue.pop_front();
			auto unzipped_data = decompressBlock(block);
			for (auto b : unzipped_data) bytes.push_back(b);
			return make_pair(block.start, block.num_alignments);
		}
		else {
			if (chromosome_trees.find(chromo) == chromosome_trees.end() ) {
				cerr << "[INFO] No data for chromosome " << chromo << endl;
				return make_pair(-1, 0);
			}
			auto tree = chromosome_trees[chromo];
			vector<RawDataInterval> overlapping;
			// find the first available coordinate
			auto first_it = tree.getFirstInterval();
			if (first_it.chromosome < 0) {
				// no intervals in a tree
				cerr << "[INFO] No data for chromosome " << chromo << endl;
				return make_pair(-1, 0);
			}
			int actual_start_coord = start_coord;
			if (start_coord < first_it.start) {
				actual_start_coord = first_it.start;
				is_transcript_start = true;
			}
			// cerr << "Load overlapping block for chr=" << chromo << ":" << actual_start_coord << endl;

			tree.findOverlapping(actual_start_coord, end_coord, overlapping);
			if (overlapping.size() > 0) {
				no_blocks = false;
				// queue up all the blocks
				for (auto b : overlapping) block_queue.push_back(b);
				// get the first block, prepare to serve it to the wrapping buffer
				RawDataInterval block = block_queue.front();
				is_transcript_start = block.isAlignedWithTranscriptStart();
				block_queue.pop_front();
				auto unzipped_data = decompressBlock(block);
				// enqueue the bytes for further use
				for (auto b : unzipped_data) bytes.push_back(b);
				return make_pair(block.start, block.num_alignments);
			}
			else {
				// no data
				no_blocks = true;
			}
		}
		return make_pair(-1, 0);
	}

	////////////////////////////////////////////////////////////////////////////
	// returns true if file is open and more bytes are available for reading
	////////////////////////////////////////////////////////////////////////////
	bool opened() {return f_in.is_open();}

	////////////////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////////////////
	bool hasMoreBytes() {
		// f_in.peek(); // peek -- will set the eof bits if reached the end of file
		// return bytes.size() > 0 || ( !no_blocks && f_in.good() ); // either have bytes in the buffer or have not reached eof
		return bytes.size() > 0 || block_queue.size() > 0;
	}

	////////////////////////////////////////////////////////////////////////////
	// stream next byte out, if no bytes available -- try to read more
	// (decompress if necessary)
	////////////////////////////////////////////////////////////////////////////
	uint8_t getNextByte() {
		if (bytes.size() < 1) {
			readMoreLZIPBlocks();
		}
		uint8_t c = bytes.front();
		bytes.pop_front();
		return c;
	}

	////////////////////////////////////////////////////////////////////////////
	vector<uint8_t> getNextNBytes(int n) {
		if (bytes.size() < n) {
			readMoreLZIPBlocks();
		}
		vector<uint8_t> local_bytes;
		while (n > 0 && bytes.size() > 0) {
			local_bytes.push_back(bytes.front());
			bytes.pop_front();
			n--;
		}
		return local_bytes;
	}

	////////////////////////////////////////////////////////////////////////////
	void popNBytes(int n) {
		if (bytes.size() < n) readMoreLZIPBlocks();

		while (n > 0 && bytes.size() > 0) {
			bytes.pop_front();
			n--;
		}
	}
};

#endif
