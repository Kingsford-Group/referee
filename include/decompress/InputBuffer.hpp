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
typedef int chromoID;




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
	int size;	// block's size (including header and trailer)
	int offset;	// offset to the block's header from the begining of the file

	MyBlock(int s, int off): size(s), offset(off) {}
};



////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
class InputBuffer {
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

		// TODO: need a bunch of interval trees -- one per chromosome, for example
		// read from the begining of the file or from the reference file

		// fill out chromosome_trees
		// RawDataInterval
		// IntervalTree tree;
		auto lzip_blocks = seek_blocks(f_in);
		int chromoID = -1;
		for (int i = 0; i < genomic_intervals->size(); i++) {
			auto interval = genomic_intervals->at(i);
			auto block = lzip_blocks[i];
			if (chromoID != interval.start.chromosome) {
				// TODO new tree
				chromoID = interval.start.chromosome;
				// tree = new tree();

				// auto chromo_span = interval.end.chromosome - interval.start.chromosome;
				// int i = 0;
				// while (i < chromo_span) {
				// 	// interval spans over chromosomes
				// 	i++;
				// }
				// tree.add(interval);
				// chromosome_trees[chromoID] = tree;
			}
			else {
				// TODO
			}
		}
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

	unordered_map<chromoID, IntervalTree<int,int> > chromosome_trees;
	// unordered_map<chromoID, IntervalTree<int,vector<uint8_t> > > unzipped_data_trees;

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
			pos -= member_size;
			blocks.emplace_back(member_size, pos);
			f_in.seekg( -member_size, f_in.cur);
		}
		reverse( blocks.begin(), blocks.end() );
		return blocks;
	}
};

#endif