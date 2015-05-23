#ifndef INPUT_STREAM_LIB_H
#define INPUT_STREAM_LIB_H

#include <vector>
#include <queue>
#include <deque>
#include <memory>

#include <fcntl.h>
#include <unistd.h>

// #include "plzip/file_index.h"

#include "tbb/concurrent_queue.h"

#include "IntervalTree.h"

using namespace std;
using namespace tbb;

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
//
//
////////////////////////////////////////////////////////////////
class InputBuffer {

	////////////////////////////////////////////////////////////////
	void createChromosomeIntervalTree(int const chromoID,
			shared_ptr<vector<TrueGenomicInterval>> genomic_intervals,
			vector<MyBlock> & lzip_blocks,
			vector<RawDataInterval> & leftovers,
			int i, int prev_i) {
		auto first_gi = genomic_intervals->at(prev_i);
		auto last_gi = genomic_intervals->at(i - 1);
		vector<RawDataInterval> chromosome_intervals;

		// add leftovers to the vector
		for (auto i : leftovers)
			chromosome_intervals.push_back(i);
		// chromosome_intervals.insert(leftovers.begin(), leftovers.end() );
		leftovers.clear();
		// 
		for (auto j = 0; j < i - prev_i; j++) {
			auto gi = genomic_intervals->at(prev_i + j);
			auto block = lzip_blocks[prev_i + j];
			// signature:
			// RawDataInterval(size_t off, size_t size, size_t ex, int chr, int s, int e): 
			if (gi.start.chromosome != gi.end.chromosome) {
				// spans several chromosomes
				auto chromo_span = gi.end.chromosome - gi.start.chromosome - 1;
				chromosome_intervals.emplace_back(
					block.offset, block.compressed_size, block.decompressed_size,
					gi.start.chromosome, gi.start.offset, chromo_max);
				chromosome_trees[gi.start.chromosome] = 
					IntervalTree<int,int>(chromosome_intervals);

				// start a new tree!
				// add default spans for all chromosomes in between
				int k = 0;
				while (k < chromo_span) {
					chromosome_intervals.clear();
					chromosome_intervals.emplace_back(
						block.offset, block.compressed_size, block.decompressed_size,
						gi.start.chromosome + k + 1, chromo_min, chromo_max);
					chromosome_trees[gi.start.chromosome + k + 1] = 
						IntervalTree<int,int>(chromosome_intervals);
					k++;
				}
				// technically, this is LEFTOVERS
				// chromosome_intervals.clear();
				leftovers.clear();
				leftovers.emplace_back(
					block.offset, block.compressed_size, block.decompressed_size,
					gi.end.chromosome, chromo_min, gi.end.offset);
				// return this interval at the end
			}
			else {
				// contained within chromosome
				chromosome_intervals.emplace_back(
					block.offset, block.compressed_size, block.decompressed_size,
					gi.start.chromosome, gi.start.offset, gi.end.offset );
			}
		}

		if (last_gi.end.chromosome == first_gi.start.chromosome) {
			chromosome_trees[chromoID] = IntervalTree<int,int>(chromosome_intervals);
			chromosome_intervals.clear();
		}
		// return chromosome_intervals;
	}

public:
	////////////////////////////////////////////////////////////////
	// default constructor: acquires file descriptor for a data stream,
	// builds a mapping between genomic coordinates and the byte blocks
	// within lzip-ed stream
	////////////////////////////////////////////////////////////////
	InputBuffer(string const & fname, shared_ptr<vector<TrueGenomicInterval>> genomic_intervals, int id, int bs = 1<<22):
		buffer_id(id),
		buffer_size(bs),
		f_in(fname.c_str(), ifstream::in)  {
			// f_in = open(fname.c_str(), ifstream::in);
			if (!f_in) {
				cerr << "[INFO] Could not open file: " << fname << endl;
			}

		// interval trees -- one per chromosome
		// fill out chromosome_trees
		auto lzip_blocks = seek_blocks(f_in);
		int chromoID = genomic_intervals->at(0).start.chromosome;
		int i, prev_i = 0;
		vector<RawDataInterval> leftovers;
		for (i = 1; i < genomic_intervals->size(); i++) {
			auto interval = genomic_intervals->at(i);
			if (chromoID != interval.start.chromosome) {
				// create new tree(s)
				// genomic_intervals[prev_i, i-1];
				createChromosomeIntervalTree(chromoID, genomic_intervals, lzip_blocks, leftovers, i, prev_i);
				// update chromo ID
				chromoID = interval.end.chromosome;
				prev_i = i;
			}
		}
		cerr << fname << " trees: " << chromosome_trees.size() << endl;
		// process remaining intervals
		// genomic_intervals[prev_i:]
		createChromosomeIntervalTree(chromoID, genomic_intervals, lzip_blocks, leftovers, i, prev_i);

	}

	////////////////////////////////////////////////////////////////
	~InputBuffer() {
		f_in.close();
	}

	////////////////////////////////////////////////////////////////
	// find all blocks overlapping this genomic interval
	// return a vector w/ block starts (in byte offsets form the beginnning of the file)
	////////////////////////////////////////////////////////////////
	deque<RawDataInterval> findOverlappingBlocks(GenomicInterval & interval) {
		deque<RawDataInterval> overlapping_blocks;

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
	void enqueueBlock(const RawDataInterval & block, concurrent_queue<shared_ptr<DataPacket>> & Q) {
		// set file cursor to the needed block (absolute position)
		f_in.seekg(block.byte_offset);

		// read bytes
		shared_ptr<vector<uint8_t>> raw_bytes(new vector<uint8_t>(block.block_size) );
		f_in.read( (char *) &(*raw_bytes)[0], block.block_size);

		shared_ptr<DataPacket> packet( new DataPacket(buffer_id, block, raw_bytes, block.decompressed_size) );
		Q.push(packet);
	}

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
		return bytes.size() > 0 || f_in.good(); // either have bytes in the buffer or have not reached eof
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
		cerr << "reading blocks from plziped stream" << endl;
		vector<MyBlock> blocks;
		// set cursor at the end of the stream
		f_in.seekg(0, f_in.end);
		// length of file
		int64_t pos = f_in.tellg();

		File_header header;
  		File_trailer trailer;
  		// scan until cursor reaches the begining of the file
  		while (f_in.tellg() > 0) {
  			// read trailer
			f_in.read( (char *) trailer.data, File_trailer::size);
			// block size
			auto member_size = trailer.member_size();
			auto data_size = trailer.data_size();
			pos -= member_size;
			blocks.emplace_back(member_size, data_size, pos);
			f_in.seekg( -member_size, f_in.cur);
		}
		reverse( blocks.begin(), blocks.end() );
		return blocks;
	}
};

#endif