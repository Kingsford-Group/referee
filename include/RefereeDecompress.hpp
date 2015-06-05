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
}

/*////////////////////////////////////////////////////////////////
Parse genomic coordinates file, return a map of vectors of intervals
grouped by type of data stream
////////////////////////////////////////////////////////////////*/
unordered_map<string,shared_ptr<vector<TrueGenomicInterval>>>
	parseGenomicIntervals(string const & fname) {
	ifstream f_in(fname);
	check_file_open(f_in, fname);
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
////////////////////////////////////////////////////////////////
GenomicInterval parseInputInterval(string const & location) {
	size_t idx = location.find(":");
	if (idx == string::npos) {
		cerr << "[ERROR] Can not parse the input interval. Expected format: chr2:5000000-100000000" << endl;
		exit(1);
	}
	int chr = stoi( location.substr(0, idx).substr(3) );
	size_t idx2 = location.find('-', idx);
	if (idx2 == string::npos) {
		cerr << "[ERROR] Can not parse the input interval. Expected format: chr2:5000000-100000000" << endl;
		exit(1);
	}
	int start_coord = stoi( location.substr(idx+1, idx2 - idx) );
	int end_coord = stoi( location.substr(idx2+1) );

	cerr << chr << " " << start_coord << " " << end_coord << endl;

	return GenomicInterval(chr, start_coord, end_coord);
}


////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
int decompressFileSequential(string const & file_name, string const & ref_file_name,
	string const & fname_out, string const & location) {

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

	if (location.size() > 0) {
		// parse location str
		GenomicInterval requested_interval = parseInputInterval(location);
		// 100mbp - 105mbp --spans blocks for has_edits
		// GenomicInterval requested_interval(0, 100000000, 105000000);
		stitchAlignmentsSerial(input_streams, requested_interval, file_name, fname_out, ref_file_name);
	}
	else {
		// keep stitching alignments while there is data available
		int read_len = 100;
		Decompressor D(file_name, fname_out, ref_file_name);
		D.decompress(read_len, input_streams, D_SEQ);
	}
}

#endif
