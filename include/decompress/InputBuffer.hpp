#ifndef INPUT_STREAM_LIB_H
#define INPUT_STREAM_LIB_H

#include <vector>
#include <queue>
#include <deque>
#include <memory>

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
//
//
////////////////////////////////////////////////////////////////
class DataPacket {
public:
	DataPacket(const int id, RawDataInterval i, shared_ptr<vector<uint8_t>> d,
		int s):
		interval(i),
		buffer_id(id),
		data(d) {}

	shared_ptr<vector<uint8_t>> getRawData() {
		return data;
	}

	int getExpectedDataSize() {return interval.decompressed_size;}

	int getBufferID() {
		return buffer_id;
	}

	RawDataInterval getInterval() {return interval;}

private:
	RawDataInterval interval;
	int buffer_id;
	// bool is_compressed = true;
	shared_ptr<vector<uint8_t>> data;
};


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


// write raw data into decoder's stream
void writeToBuffer(LZ_Decoder * const decoder, shared_ptr<vector<uint8_t>> raw_data) {
	
}

////////////////////////////////////////////////////////////////
//
// http://www.nongnu.org/lzip/manual/lzlib_manual.html#Examples
// example 8
//
////////////////////////////////////////////////////////////////
vector<uint8_t> unzipData(shared_ptr<vector<uint8_t>> raw_data, int const new_data_size) {
	cerr << "Unzipping LZIP block expected size: " << new_data_size << endl;
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
	// cerr << "Can write " << chunk_size << " need to write " << need_to_write << endl;


	int nz = 0;
	// cerr << "magic: " << (*raw_data)[0] << (*raw_data)[1] << (*raw_data)[2] << 
	// 	(*raw_data)[3] << 
	// 	" version: " << (int)(*raw_data)[4] << 
	// 	" coded_dict_size: " << (int)(*raw_data)[5] << endl;
	// for (auto c : *raw_data)
	// 	nz += (c > 0);
	// cerr << "block has: " << nz << " non-zero entries" << endl;

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
		cerr << "wrote " << written << " ";

		int read = LZ_decompress_read( decoder, (uint8_t *) &unzipped_data[bytes_read], new_data_size - bytes_read );
		bytes_read += read;
		cerr << "read " << read << " ";

		chunk_size = LZ_decompress_write_size( decoder );
	}
	cerr << "total bytes written: " << wrote << " total bytes read: " << bytes_read << endl;
	assert( wrote >= need_to_write );
	assert(bytes_read == new_data_size);

	// int read = LZ_decompress_read( decoder, (uint8_t *) &unzipped_data[bytes_read], new_data_size - bytes_read );
	// cerr << " read " << read << " ";
	// if( read < 0 ) {
	// 	if( LZ_decompress_errno( decoder ) == LZ_header_error ) {
	// 		cerr << "sync to memebr ";
	// 		LZ_decompress_sync_to_member( decoder );
	// 		int read = LZ_decompress_read( decoder, (uint8_t *) &unzipped_data[bytes_read], new_data_size - bytes_read );
	// 		cerr << " read2 " << read << " ";
 //        	// trailing_garbage_found = true;
 //        }
 //        else {
	// 		cerr << "[ERROR] LZIP decompress error" << endl;
	// 		exit(1);
	// 	}
	// }
	
	LZ_decompress_finish( decoder );

	auto retval = LZ_decompress_finished( decoder );
	// destroys all internal data
	LZ_decompress_close( decoder );

	nz = 0;
	for (auto c : unzipped_data) nz += (c > 0);

	cerr << "Data: nz=" << nz << " ";
	for (int i = 0; i < 200; i++) cerr << unzipped_data[i];
		cerr << endl;

	// TODO: might be better to keep around a decoder if set up is expensive
	// LZ_decompress_reset( decoder );	// prepare for new member

	return unzipped_data;
}


////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
class InputBuffer {

	void createTree(int const chromo, 
		vector<RawDataInterval> const & chromo_intervals, 
		unordered_map<chromo_id_t, IntervalTree<int,int> > & chromosome_trees) {
		chromosome_trees[chromo] = 
				IntervalTree<int,int>(chromo_intervals);
		cerr << "Chromosome tree: " << chromosome_trees.size();
	}

	void createChromosomeIntervalTree2(
		shared_ptr<vector<TrueGenomicInterval>> genomic_intervals, 
		vector<MyBlock> & lzip_blocks, 
		unordered_map<chromo_id_t, IntervalTree<int,int> > & chromosome_trees) {

		// int range_start = 0, range_end = 0;
		int prev_chromo = genomic_intervals->at(0).start.chromosome;
		vector<RawDataInterval> chromo_intervals;
		for (int i = 0; i < genomic_intervals->size(); i++) {
			auto interval = genomic_intervals->at(i);
			if (prev_chromo != interval.start.chromosome) {
				// create a tree for intervals in [range_start, range_end]
				createTree(prev_chromo, chromo_intervals, chromosome_trees);
				// start a new range at this index
				chromo_intervals.clear();
				// range_start = i;
				prev_chromo = interval.start.chromosome;
			}
			
			// interval.start.chromosome is prev_chromo
			if (interval.start.chromosome == interval.end.chromosome) {
				// start and end chromosome are the same and equal prev_chromo
				// range_end = i;
				chromo_intervals.emplace_back(
					block.offset, block.compressed_size, block.decompressed_size,
					gi.start.chromosome, gi.start.offset, gi.end.offset);
				// moving on...
			}
			else {
				// figure out how many chromos this interval spans
				auto chromo_span = interval.end.chromosome - interval.start.chromosome - 1;
				// split the interval as many times as needed (at least 2 pieces)
				// first create a tree for chromo_intervals w/ a starting piece of the current interval
				chromo_intervals.emplace_back(
					block.offset, block.compressed_size, block.decompressed_size,
					gi.start.chromosome, gi.start.offset, chromo_max);
				// then create tree(s) for every chromo that is in between the start and the end
				chromo_intervals.clear();
				int k = 0;
				while (k < chromo_span) {
					chromo_intervals.emplace_back(
						block.offset, block.compressed_size, block.decompressed_size,
						gi.start.chromosome + k + 1, chromo_min, chromo_max);
					createTree(gi.start.chromosome + k + 1, chromo_intervals, chromosome_trees);
					chromo_intervals.clear();
					k++;
				}
				// then add the leftover piece to the chromo_intervals
				chromo_intervals.emplace_back(
					block.offset, block.compressed_size, block.decompressed_size,
					gi.end.chromosome, chromo_min, gi.end.offset);
				prev_chromo = gi.end.chromosome;
			}
			
		}
		// create a tree for intervals in chromo_intervals
		createTree(prev_chromo, chromo_intervals, chromosome_trees);
	}


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
			// f_in = open(fname.c_str(), ifstream::in);
			if (!f_in) {
				cerr << "[INFO] Could not open file: " << fname << endl;
				return;
			}

		// interval trees -- one per chromosome
		// fill out chromosome_trees
		auto lzip_blocks = seek_blocks(f_in);
		cerr << fname << " blocks: " << lzip_blocks.size() << endl;
		cerr << "BUILDING TREES" << endl;

		createChromosomeIntervalTree2(genomic_intervals, lzip_blocks, chromosome_trees);
		cerr << " Total trees: " << chromosome_trees.size() << endl;
	}

	////////////////////////////////////////////////////////////////
	~InputBuffer() {
		f_in.close();
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	int loadOverlappingBlock(int chromo, int start_coord,
		bool & is_transcript_start) {
		
		bytes.clear();
		auto tree = chromosome_trees[chromo];
		vector<RawDataInterval> overlapping;

		// todo: find the first available coordinate
		auto first_it = tree.getFirstInterval();
		if (first_it.chromosome < 0) {
			// no intervals in a tree
			return -1;
		}
		if (start_coord < first_it.start) {
			start_coord = first_it.start;
			is_transcript_start = true;
		}

		cerr << "Load overlapping block for chr=" << chromo << ":" << start_coord << endl;
		tree.findOverlapping(start_coord, start_coord + 1, overlapping);
		cerr << "Num overlapping blocks: " << overlapping.size() << endl;
		if (overlapping.size() > 0) {
			no_blocks = false;
			// decompress the first block, prepare to serve it to the wrapping buffer
			RawDataInterval block = overlapping[0];
			// read block bytes from the LZ stream
			f_in.seekg(0, f_in.beg);
			f_in.seekg(block.byte_offset);
			shared_ptr<vector<uint8_t>> raw_bytes(new vector<uint8_t>(block.block_size) );
			f_in.read( (char *) &(*raw_bytes)[0], block.block_size);
			cerr << "Read " << block.block_size << " bytes" << endl;
			// for (auto c : *raw_bytes) cerr << (int)c;
			// cerr << "------";
			auto unzipped_data = unzipData(raw_bytes, block.decompressed_size);
			cerr << "Expected: " << block.decompressed_size << " bytes; got: " << unzipped_data.size() << endl;
			for (auto b : unzipped_data)
				bytes.push_back(b);
			cerr << block.start << endl;
			return block.start;
		}
		else {
			// no data
			no_blocks = true;
		}
		return start_coord;
	}

	////////////////////////////////////////////////////////////////
	// find all blocks overlapping this genomic interval
	// return a vector w/ block starts (in byte offsets form the beginnning of the file)
	////////////////////////////////////////////////////////////////
	deque<RawDataInterval> findOverlappingBlocks(GenomicInterval & interval) {
		deque<RawDataInterval> overlapping_blocks;

		// cerr << "Tree cnt: " << chromosome_trees.size() << endl;
		int chromo = interval.chromosome;
		auto tree = chromosome_trees[chromo];
		vector<RawDataInterval> overlapping_intervals;
		// TODO: overload method to allow deque
		tree.findOverlapping(interval.start, interval.stop, overlapping_intervals);

		for (auto intvl : overlapping_intervals) {
			// create a block w/ byte offset, chromosome, appropriate start, end
			overlapping_blocks.emplace_back(intvl);
				// intvl.value, block_size,
				// 		decompressed_size, chromo, intvl.start, intvl.stop);
		}
		return overlapping_blocks;
	}

	//////////////////////////////////////////////////////////////// 
	// read bytes for this block from the file
	// create a DataPacket, enqueue it
	// should only overlap a single block
	////////////////////////////////////////////////////////////////
	// void enqueueBlock(const RawDataInterval & block, shared_ptr<concurrent_queue<shared_ptr<DataPacket>>> Q) {
	// 	cerr << name.substr(name.size() - 10) << ": chr" << block.chromosome << ":" <<
	// 		block.start << "-" << block.stop << endl;
	// 	// set file cursor to the needed block (absolute position)
	// 	f_in.seekg(block.byte_offset);

	// 	// read bytes
	// 	shared_ptr<vector<uint8_t>> raw_bytes(new vector<uint8_t>(block.block_size) );
	// 	f_in.read( (char *) &(*raw_bytes)[0], block.block_size);

	// 	shared_ptr<DataPacket> packet( new DataPacket(buffer_id, block, raw_bytes, block.decompressed_size) );
	// 	Q->push(packet);
	// }

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	// void addData(RawDataInterval & interval, vector<uint8_t> & unzipped_data) {
		// auto tree = unzipped_data_trees[interval.chromosome];
		// tree.add(interval, unzipped_data);
	// }

	////////////////////////////////////////////////////////////////
	// returns true if file is open and more bytes are available for reading
	////////////////////////////////////////////////////////////////
	bool opened() {return f_in.is_open();}


	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	bool hasMoreBytes() {
		f_in.peek(); // peek -- will set the eof bits if reached the end of file
		return bytes.size() > 0 || ( !no_blocks && f_in.good() ); // either have bytes in the buffer or have not reached eof
	}

	////////////////////////////////////////////////////////////////
	// make this a blocking function -- if stuff not yet available
	// wait for LZIP to decompress it
	// potential for deadlock?
	////////////////////////////////////////////////////////////////
	uint8_t getNextByte() {
		if (bytes.size() < 1) {
			readMore();
		}
		uint8_t c = bytes.front();
		bytes.pop_front();
		// current_index++;
		return c;
	}

	////////////////////////////////////////////////////////////////
	vector<uint8_t> getNextNBytes(int n) {
		if (bytes.size() < n) {
			readMore();
		}
		vector<uint8_t> local_bytes;
		// for (int i = 0; i < n && i < bytes.size(); i++) {
		while (n > 0 && bytes.size() > 0) {
			local_bytes.push_back(bytes.front());
			bytes.pop_front();
			// current_index++;
			n--;
		}
		return local_bytes;
	}

	void popNBytes(int n) {
		if (bytes.size() < n) readMore();

		while (n > 0 && bytes.size() > 0) {
			bytes.pop_front();
			// current_index++;
			n--;
		}
	}

	// int current_index = 0;
////////////////////////////////////////////////////////////////
private:

	bool no_blocks = false;

	// file stream name
	string name;

	unordered_map<chromo_id_t, IntervalTree<int,int> > chromosome_trees;
	// unordered_map<chromo_id_t, IntervalTree<int,vector<uint8_t> > > unzipped_data_trees;

	int buffer_id;

	ifstream f_in;

	deque<uint8_t> bytes;

	// File_index file_index;

	int buffer_size;

	////////////////////////////////////////////////////////////////
	void readMore() {
		// TODO: when need to read more have a worker thread read more bytes from an
		// appropriate stream
		// decompress 
		// add decompressed bytes to this buffer
		vector<uint8_t> chunk(buffer_size, 0);
		f_in.read( (char *) &chunk[0], buffer_size);
		int actually_read = f_in.gcount();
		if (actually_read < buffer_size) {
			// cerr << "Requested " << buffer_size << ", got " << actually_read << endl;
		}
		for (auto c : chunk)
			bytes.push_back(c);
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
		// shift by trailer::size
		f_in.seekg( -File_trailer::size, f_in.cur);
		// total += File_trailer::size;
		cerr << "file len: " << file_len << " ";

		File_header header;
  		File_trailer trailer;
  		// scan until cursor reaches the begining of the file
  		while (pos > 0) {
  			// read trailer
			f_in.read( (char *) trailer.data, File_trailer::size);
			// block size
			auto member_size = trailer.member_size();
			auto data_size = trailer.data_size();
			// cerr << "m: " << member_size << " d: " << data_size << " ";
			pos -= member_size;
			total += member_size;
			blocks.emplace_back(member_size, data_size, pos);
			f_in.seekg( -member_size - File_trailer::size, f_in.cur);
		}
		cerr << "..." << name.substr(name.size() - 10) << 
			" File len: " << file_len << " byte seen: " << total << endl;
		assert(total == file_len);
		reverse( blocks.begin(), blocks.end() );
		return blocks;
	}
};

#endif
