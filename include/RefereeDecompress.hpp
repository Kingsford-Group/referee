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

#include "decompress/Decompressor.hpp"

////////////////////////////////////////////////////////////////
//
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
	uint8_t options = 1;
	D.decompressInterval(requested_interval, read_len, input_streams, options);

	cerr << "Covered the interval (or no more alignments within the interval)" << endl;
}

/*////////////////////////////////////////////////////////////////
Parse genomic coordinates file, return a map of vectors of intervals
grouped by type of data stream
////////////////////////////////////////////////////////////////*/
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
		auto second_space = line.find(' ', space + 1);
		unsigned long num_alignments = stoul(line.substr(space + 1, second_space));
		map[suffix]->emplace_back( line.substr(second_space + 1), num_alignments );
	}
	f_in.close();
	return map;
}


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
	// TODO: don't need the buffer_map anymore
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
	GenomicInterval requested_interval(0, 0, 10000000);
	// 100mbp - 105mbp --spans blocks for has_edits
	// GenomicInterval requested_interval(0, 100000000, 105000000);

	stitchAlignmentsSerial(input_streams, requested_interval, file_name, fname_out, ref_file_name);
}

#endif