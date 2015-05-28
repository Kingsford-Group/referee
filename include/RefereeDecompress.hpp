#ifndef REFEREE_DECOMPRESS_LIB_H
#define  REFEREE_DECOMPRESS_LIB_H

#include <fcntl.h>
#include <unistd.h>
#include <lzip.h>
#include <compress.h>
#include <chrono>
#include <thread>
#include <sys/types.h>


#include "file_index.h"

// #include "tbb/concurrent_queue.h"
// #include "tbb/concurrent_vector.h"
// #include "tbb/concurrent_unordered_map.h"



#include "decompress/Decompressor.hpp"

// using namespace tbb;

////////////////////////////////////////////////////////////////
//
// find blocks overlapping the interval
// add them to a queue to be decompressed checking that the queue does not
// exceed a certain size limit
//
////////////////////////////////////////////////////////////////
// void enqueueBlocks(
// 	shared_ptr<concurrent_queue<shared_ptr<DataPacket>>> Q,
//         shared_ptr<vector<GenomicInterval>> intervals,
//         shared_ptr<vector<shared_ptr<InputBuffer>>> inputs,
//         shared_ptr<bool> done) {
// 	// do not need to sync since using a lock-free Q
// 	for (auto interval : *intervals) {
// 		// cerr << interval.chromosome << ":" << interval.start << "-" << interval.stop << endl;
// 		// find blocks in all relevant input streams that overlap this interval
// 		bool have_blocks = false;

// 		vector<deque<RawDataInterval>> block_starts;
// 		// cerr << "Num of inputs: " << inputs.size() << endl;
// 		for (auto input_buffer : *inputs) {
// 			// TODO: does push use move or deep copy?
// 			block_starts.push_back( input_buffer->findOverlappingBlocks(interval) );
// 			have_blocks = have_blocks || (block_starts.back().size() > 0);
// 		}
// 		cerr << "Have blocks: " << have_blocks << endl;

// 		// keep adding compressed blocks of data to the queue in order of their appearance 
// 		// in their respecting plzip streams
// 		while (have_blocks) {
// 			have_blocks = false;
// 			for (int i = 0; i < inputs->size(); i++) {
// 				have_blocks = have_blocks || (block_starts[i].size() > 0);
// 				if (block_starts[i].size() > 0) {
// 					auto block = block_starts[i].front(); // get the start byte for the next block
// 					block_starts[i].pop_front(); // pop it off the front
// 					// cerr << "enqueuing ";
// 					inputs->at(i)->enqueueBlock(block, Q);
// 				}
// 			}
// 		}
// 	}

// 	// #pragma omp critical (CRITICAL)
//     // {
//     	*done = true;
//     	cerr << "Done enqueuing all the needed blocks" << endl;
//     // }
// }



////////////////////////////////////////////////////////////////
//
// pull things off the queue until it is empty
// decompress them using lzip functions
// add decompressed data to the appropriate buffer
//
////////////////////////////////////////////////////////////////
// void lzipDecompress(
// 	shared_ptr<concurrent_queue<shared_ptr<DataPacket>>> Q,
//         shared_ptr<bool> done) {

// 	shared_ptr<DataPacket> data_packet;
// 	while (!done || !Q->empty()) {
//         if (Q->try_pop(data_packet) ) {
// 		cerr << "popped ";
//         	// once data is decompressed -- hand it back to packet and let
//         	// packet add it to its InputBuffer
//         	auto unzipped_data = unzipData(data_packet->getRawData(), data_packet->getExpectedDataSize() );
//         	auto buffer_id = data_packet->getBufferID();
//         	auto interval = data_packet->getInterval();
// 			// buffer_map[buffer_id]->addData(interval, unzipped_data);
// 			interval.setDecompressedData(unzipped_data);
// 		}
// 	}

// 	//#pragma omp critical (CRITICAL)
// 	//{
//     	// cerr << "Done decompressing all the needed blocks" << endl;
// 	//}
// }


////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
// void stitchAlignments(shared_ptr<InputStreams> input_streams,
// 		shared_ptr<vector<GenomicInterval>> intervals,
// 		shared_ptr<string> output_name, shared_ptr<string> ref_name,
// 	        shared_ptr<bool> done) {


// 	// keep stitching alignments while there is data available
// 	Decompressor D(*output_name, *ref_name);
// 	// will wait for InputBuffers to get populated
// 	// and also -- wait for "done" to be true
// 	// D.decompress(input_streams, done);

// 	for (auto interval : *intervals) {
// 		// TODO: request data for this interval from the input_streams
// 		// cerr << "TODO";
// 	}

// 	// #pragma omp critical (CRITICAL)
// 	// {
// 		cerr << "done stitching alignments" << endl;
// 	// }
// }


////////////////////////////////////////////////////////////////
void stitchAlignmentsSerial(
	InputStreams & input_streams,
	GenomicInterval & requested_interval, 
	string const & input_fname, 
	string const & output_name, 
	string const & ref_name) {
	// keep stitching alignments while there is data available
	int read_len = 100;
	Decompressor D(input_fname, output_name, ref_name);
	D.decompressInterval(requested_interval, read_len, input_streams);

	cerr << "Covered the interval (or no more alignments within the interval)" << endl;
}

/*
Parse genomic coordinates file, return a map of vectors of intervals
grouped by type of data stream
*/
unordered_map<string,shared_ptr<vector<TrueGenomicInterval>>> 
	parseGenomicIntervals(string const & fname) {
		ifstream f_in(fname);
		if (!f_in) {
			cerr << "[ERROR] Could not find a mapping to genomic coordinates" << endl;
			exit(1);
		}
		string line;
		unordered_map<string, shared_ptr<vector<TrueGenomicInterval>> > map;
		while ( getline(f_in, line) ) {
			auto space = line.find(' ');
			string suffix = line.substr(0, space);
			if (map.find(suffix) == map.end()) {
				shared_ptr<vector<TrueGenomicInterval>> intervals(new vector<TrueGenomicInterval>());
				map[suffix] = intervals;
			}
			map[suffix]->emplace_back( line.substr(space+1) );
			// map[suffix]->back().print();
		}
		f_in.close();
		return map;
	}


////////////////////////////////////////////////////////////////
//
// Read file indices, assign threads to parse blocks, wait for
// threads to finish (join) and clean up.
//
////////////////////////////////////////////////////////////////
// int decompressFile(string const & file_name, string const & ref_file_name, 
// 	string const & fname_out, const int num_workers) {

// 	// set up inputs
// 	unordered_map<string,shared_ptr<vector<TrueGenomicInterval>>> all_intervals = 
// 		parseGenomicIntervals("genomic_intervals.txt");
// 	InputStreams input_streams;

// 	int buffer_size = pow(2, 24); // 16Mb
// 	int buffer_id = 0;
// 	shared_ptr<unordered_map<int, shared_ptr<InputBuffer>>> buffer_map(
// 		new unordered_map<int, shared_ptr<InputBuffer>>());

// 	shared_ptr<vector<shared_ptr<InputBuffer>>> input_buffers(
// 		new vector<shared_ptr<InputBuffer>>() );
// 	for (auto pair : all_intervals) {
// 		auto suffix = pair.first;
// 		auto intervals = pair.second;
// 		// cerr << suffix << endl;
// 		shared_ptr<InputBuffer> ib(new InputBuffer(file_name + suffix,
// 			intervals, buffer_size, buffer_id));
// 		if ( !ib->opened() ) {
// 			ib.reset();
// 			continue;
// 		}
// 		input_buffers->push_back(ib);
// 		buffer_map->emplace(buffer_id++, ib);

// 		if (suffix.compare(".offs.lz") == 0) {
// 			input_streams.offs = shared_ptr<OffsetsStream>(new OffsetsStream(ib) );
// 		}
// 		//else if (suffix.compare(".edits.lz") == 0 ) {
// 		//	input_streams.edits = shared_ptr<EditsStream>(new EditsStream(ib) );
// 		//}
// 		else if (suffix.compare("left_clip.lz") == 0) {
// 			input_streams.left_clips = shared_ptr<ClipStream>(new ClipStream(ib) );
// 		}
// 		else if (suffix.compare("right_clip.lz") == 0) {
// 			input_streams.right_clips = shared_ptr<ClipStream>(new ClipStream(ib) );
// 		}
// 		// TODO: add more inputs
// 	}
// 	cerr << "Initiated streams" << endl;

// 	// input intervals: the ones that need to be decompressed
// 	shared_ptr<vector<GenomicInterval>> intervals( new vector<GenomicInterval>() );
// 	GenomicInterval test_interval(0, 0, 1000000);
// 	intervals->push_back(test_interval);

// 	// flag indicating whether there are more packets coming or not
// 	bool done = false;

// 	// set up multiple OpenMP threads
// 	//omp_set_dynamic(true);
// 	//auto max_threads = omp_get_num_procs();
// 	//omp_set_num_threads( (int)min( (int)num_workers, max_threads) );
// 	cerr << "Threads:\t" << (int)num_workers << endl;

// 	// launch enqueueing thread
// 	// everything has to be a pointer
// 	thread t1(foo, intervals);
// 	std::thread q_thread(enqueueBlocks, Q, intervals, input_buffers, make_shared<bool>(done) );

// 	// launch thread to assembly alignments
// 	std::thread s_thread(stitchAlignments, make_shared<InputStreams>(input_streams), intervals,
// 		make_shared<string>(fname_out), make_shared<string>(ref_file_name), make_shared<bool>(done) );
// 	// launch other threads to decompress blocks
// 	vector<thread> worker_threads(num_workers-2);

// 	for (int i = 0; i < num_workers-2; i++) {
// 		worker_threads[i] = thread(lzipDecompress, Q, make_shared<bool>(done) );
// 	}

// 	// join threads
// 	q_thread.join();
// 	for (int i = 0; i < num_workers-2; i++) worker_threads[i].join();
// 	s_thread.join();

// 	// all threads finish and join here
// 	cerr << "Done, done, done!" << endl;
// }



////////////////////////////////////////////////////////////////
//
// 
//
////////////////////////////////////////////////////////////////
int decompressFileSequential(string const & file_name, string const & ref_file_name, 
	string const & fname_out) {

	// set up inputs
	unordered_map<string,shared_ptr<vector<TrueGenomicInterval>>> all_intervals = 
		parseGenomicIntervals("genomic_intervals.txt");
	InputStreams input_streams;

	int buffer_size = pow(2, 24); // 16Mb
	int buffer_id = 0;
	unordered_map<int, shared_ptr<InputBuffer>> buffer_map;

	//////////////////////////////////////////////////////////
	// initialize input buffers
	//////////////////////////////////////////////////////////
	// vector<shared_ptr<InputBuffer>> input_buffers;
	for (auto pair : all_intervals) {
		auto suffix = pair.first;
		auto intervals = pair.second;
		if (suffix.compare(".offs.lz") == 0) {
			shared_ptr<InputBuffer> offset_ib(new InputBuffer(file_name + suffix, 
				intervals, buffer_size, buffer_id));
			// input_buffers.push_back(offset_ib);
			buffer_map.emplace(buffer_id++, offset_ib);

			input_streams.offs = shared_ptr<OffsetsStream>(new OffsetsStream(offset_ib) );
		}
		else if (suffix.compare(".edits.lz") == 0) {
			shared_ptr<InputBuffer> edits_ib(new InputBuffer(file_name + suffix, 
				intervals, buffer_size, buffer_id));
			buffer_map.emplace(buffer_id++, edits_ib);

			shared_ptr<InputBuffer> has_edits_ib(new InputBuffer(file_name + ".has_edits.lz", 
				all_intervals[".has_edits.lz"], buffer_size, buffer_id));
			buffer_map.emplace(buffer_id++, has_edits_ib);

			input_streams.edits = shared_ptr<EditsStream>(new EditsStream(edits_ib, has_edits_ib) );
		}
	}

	// input interval: the one that needs to be decompressed
	GenomicInterval requested_interval(0, 0, 1000000);

	stitchAlignmentsSerial(input_streams, requested_interval, file_name, fname_out, ref_file_name);
}

#endif