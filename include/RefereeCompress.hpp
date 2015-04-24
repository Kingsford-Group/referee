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
	int num_parsing_threads;
	Packet_courier * courier;
	bool aligned_seq_only;
	bool unique_seq_only;
};

////////////////////////////////////////////////////////////////
void * parseSAM( void * pa) {
	// cerr << "parser" << endl;
	const Parser_args & tmp = *(Parser_args *)pa;

	Packet_courier * courier = tmp.courier;
	Output_args outs = tmp.output;
	// will write out a BAM/SAM header
	Compressor c(tmp.file_name, tmp.num_parsing_threads, outs, tmp.aligned_seq_only, tmp.unique_seq_only);
	if (c.failed() ) {
		cerr << "[INFO] Terminating. " << endl;
		courier->finish();
	}
	else {
		auto start_time = chrono::system_clock::now();
		c.compress();
		auto end_time = chrono::system_clock::now();
		cerr << "Parsing wall time: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << "us" << endl;
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
Output_args initializeOutputStreams(string & name_prefix, bool aligned_seq_only, bool unique_seq_only, Packet_courier * courier) {
	Output_args oa;
	oa.offsets_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (name_prefix + ".offs.lz").c_str(), 1<<22, 20) );
	oa.edits_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (name_prefix + ".edits.lz").c_str() ) );
	oa.left_clips_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (name_prefix + ".left_clip.lz").c_str(), 1<<22, 20) );
	oa.right_clips_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (name_prefix + ".right_clip.lz").c_str(), 1<<22, 20) );
	if (!aligned_seq_only && !unique_seq_only) {
		oa.flags_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (name_prefix + ".flags.lz").c_str() ) );
		oa.ids_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (name_prefix + ".ids.lz").c_str(), 3 << 20,  12 ) );
		oa.opt_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (name_prefix + ".opt.lz").c_str() ) );
		// oa.quals_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (name_prefix + ".quals.lz").c_str() ) );
		oa.quals_buf = shared_ptr<QualityCompressor>(new QualityCompressor(courier, name_prefix.c_str(), 0.02, 200000, 4 ) );
	}
	oa.unaligned_buf = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (name_prefix + ".unaligned.lz").c_str() ) );
	
	return oa;
};

////////////////////////////////////////////////////////////////
void compressFile(string & file_name, const int num_workers, bool aligned_seq_only, bool unique_seq_only) {
	int dictionary_size = 1<<23;
	int match_len_limit = 36; // equivalent to -6 option

	// const int num_workers = num_threads - 1;//max(2, num_threads - 1); // one for parsing
	const int slots_per_worker = 20;
	const int num_slots =
		( ( num_workers > 1 ) ? num_workers * slots_per_worker : 1 );

    Packet_courier courier(num_workers, num_slots);

	// open output streams
	Output_args output_args = initializeOutputStreams(file_name, aligned_seq_only, unique_seq_only, &courier);

	// from plzip library implementation
	// initialize parsing thread -- pass courier, FDs for output
	Parser_args parser_args;
	parser_args.output = output_args;
	parser_args.file_name = file_name;
	parser_args.courier = &courier;
	parser_args.aligned_seq_only = aligned_seq_only;
	parser_args.unique_seq_only = unique_seq_only;

	pthread_t * parser_thread = new pthread_t();
	int errcode = pthread_create( parser_thread, 0, parseSAM, &parser_args );
	if( errcode ) { 
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
