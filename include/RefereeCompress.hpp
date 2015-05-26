#ifndef REFEREE_COMPRESS_LIB_H
#define  REFEREE_COMPRESS_LIB_H

#include <fcntl.h>
#include <unistd.h>
#include <lzip.h>
#include <compress.h>
#include <chrono>

#include "compress/Compressor.hpp"

struct Parser_args {
	Output_args output;
	string file_name;
	string ref_file_name;
	int num_parsing_threads;
	Packet_courier * courier;
	bool seq_only;
	bool discard_secondary_alignments;
};

////////////////////////////////////////////////////////////////
void * parseSAM( void * pa) {
	// cerr << "parser" << endl;
	const Parser_args & tmp = *(Parser_args *)pa;

	Packet_courier * courier = tmp.courier;
	Output_args outs = tmp.output;
	// will write out a BAM/SAM header
	Compressor c(tmp.file_name, tmp.ref_file_name, tmp.num_parsing_threads, outs, 
			tmp.seq_only, tmp.discard_secondary_alignments);
	if (c.failed() ) {
		cerr << "[INFO] Terminating. " << endl;
		courier->finish();
	}
	else {
		auto start_time = chrono::system_clock::now();
		c.compress();
		// finished parsing SAM -- might have data remaining in the buffers
		auto end_time = chrono::system_clock::now();
		cerr << "Parsing wall time: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << " us" << endl;
		// flush all output buffers (get rid of remaining packets)
		outs.flush();
		auto flush_time = chrono::system_clock::now();
		cerr << "Flushing wall time: " << chrono::duration_cast<chrono::microseconds>(flush_time - end_time).count() << "us" << endl;
		// let courier know that no more packages are coming
		courier->finish();
		// cerr << "finished parsing" << endl;
	}
};

////////////////////////////////////////////////////////////////
Output_args initializeOutputStreams(string & name_prefix, bool seq_only, 
		bool discard_secondary_alignments, Packet_courier * courier) {
	shared_ptr<ofstream> intervals(new ofstream("genomic_intervals.txt"));
	Output_args oa(seq_only);
	oa.offsets_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, intervals, name_prefix, ".offs.lz", 1<<22, 20) );
	oa.edits_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, intervals, name_prefix, ".edits.lz" ) );
	oa.left_clips_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, intervals, name_prefix, ".left_clip.lz", 1<<22, 20) );
	oa.right_clips_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, intervals, name_prefix, ".right_clip.lz", 1<<22, 20) );
	oa.unaligned_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, name_prefix, ".unaligned.lz" ) );
	if (!seq_only) {
		oa.flags_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, intervals, name_prefix, ".flags.lz" ) );
		oa.ids_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, intervals, name_prefix, ".ids.lz", 3 << 20,  12 ) );
		oa.opt_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, intervals, name_prefix, ".opt.lz" ) );
		// oa.quals_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, intervals, (name_prefix + ".quals.lz").c_str() ) );
		oa.quals_buf = shared_ptr<QualityCompressor>(new QualityCompressor(courier, intervals, name_prefix.c_str(), 0.05, 200000, 4 ) );
	}
	return oa;
};

////////////////////////////////////////////////////////////////
void compressFile(string & file_name, string const & ref_file_name, const int num_workers, 
	bool seq_only, bool discard_secondary_alignments) {
	int dictionary_size = 1<<23;
	int match_len_limit = 36; // equivalent to -6 option

	// const int num_workers = num_threads - 1;//max(2, num_threads - 1); // one for parsing
	const int slots_per_worker = 20;
	const int num_slots =
		( ( num_workers > 1 ) ? num_workers * slots_per_worker : 1 );

	Packet_courier courier(num_workers, num_slots);

	// open output streams
	Output_args output_args = initializeOutputStreams(file_name, seq_only, discard_secondary_alignments, &courier);

	// from plzip library implementation
	// initialize parsing thread -- pass courier, FDs for output
	Parser_args parser_args;
	parser_args.output = output_args;
	parser_args.file_name = file_name;
	parser_args.ref_file_name = ref_file_name;
	parser_args.courier = &courier;
	parser_args.seq_only = seq_only;
	parser_args.discard_secondary_alignments = discard_secondary_alignments;

	pthread_t * parser_thread = new pthread_t();
	int errcode = pthread_create( parser_thread, 0, parseSAM, &parser_args );
	if ( errcode ) { 
		show_error( "Can't create parser thread", errcode ); cleanup_and_fail(); 
	}

	// initialize worker threads
	Worker_arg worker_arg;
	worker_arg.courier = &courier;
	worker_arg.dictionary_size = dictionary_size;
	worker_arg.match_len_limit = match_len_limit;

	pthread_t * worker_threads = new( std::nothrow ) pthread_t[num_workers];
	if( !worker_threads ) { 
		cleanup_and_fail(); 
	}
	for( int i = 0; i < num_workers; ++i ) {
		errcode = pthread_create( worker_threads + i, 0, cworker, &worker_arg );
		if( errcode ) { 
			show_error( "Can't create worker threads", errcode ); cleanup_and_fail(); 
		}
	}

	// concurrently wait for threads to return compressed packets; write them to disk
	muxer(courier);

	// join worker threads
	for( int i = num_workers - 1; i >= 0; --i ) {
		errcode = pthread_join( worker_threads[i], 0 );
		if( errcode ) { 
			show_error( "Can't join worker threads", errcode ); cleanup_and_fail(); 
		}
		else {
			// cerr << "thread " << i << " joined" << endl;
		}
	}
	delete[] worker_threads;

	// join the parser thread
	errcode = pthread_join( *parser_thread, 0 );
	if( errcode ) {
		show_error( "Can't join parser thread", errcode ); cleanup_and_fail(); 
	}

	// destructor in OutputBuffer will close file output streams
};

#endif
