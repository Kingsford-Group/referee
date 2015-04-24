/*  Plzip - Parallel compressor compatible with lzip
    Copyright (C) 2009 Laszlo Ersek.
    Copyright (C) 2009, 2010, 2011, 2012, 2013 Antonio Diaz Diaz.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define _FILE_OFFSET_BITS 64

#include <algorithm>
#include <cerrno>
#include <climits>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <string>
#include <vector>
#include <iostream>
#include <pthread.h>
#include <stdint.h>
#include <unistd.h>
#include <lzlib.h>

#include "compress.h"

#ifndef LLONG_MAX
#define LLONG_MAX  0x7FFFFFFFFFFFFFFFLL
#endif


// terminate the program
void cleanup_and_fail( const int retval) {
  exit(retval);
};
void show_error( const char * const msg, const int errcode,
                 const bool help ) {
  std::cerr << "[ERROR " << errcode << "]: " << msg << std::endl;
                 }
void internal_error( const char * const msg ) {
  std::cerr << "[INFO] " << msg << std::endl;
}
void show_progress( const int packet_size,
                    const Pretty_print * const p,
                    const struct stat * const in_statsp ) {
  std::cerr << "[REDACTED]" << std::endl;
}


// Returns the number of bytes really read.
// If (returned value < size) and (errno == 0), means EOF was reached.
//
int readblock( const int fd, uint8_t * const buf, const int size )
  {
  int rest = size;
  errno = 0;
  while( rest > 0 )
    {
    const int n = read( fd, buf + size - rest, rest );
    if( n > 0 ) rest -= n;
    else if( n == 0 ) break;				// EOF
    else if( errno != EINTR && errno != EAGAIN ) break;
    errno = 0;
    }
  return size - rest;
  }


// Returns the number of bytes really written.
// If (returned value < size), it is always an error.
//
int writeblock( const int fd, const uint8_t * const buf, const int size )
  {
  int rest = size;
  errno = 0;
  while( rest > 0 )
    {
    const int n = write( fd, buf + size - rest, rest );
    if( n > 0 ) rest -= n;
    else if( n < 0 && errno != EINTR && errno != EAGAIN ) break;
    errno = 0;
    }
  return size - rest;
  }


void xinit( pthread_mutex_t * const mutex )
  {
  const int errcode = pthread_mutex_init( mutex, 0 );
  if( errcode )
    { show_error( "pthread_mutex_init", errcode ); cleanup_and_fail(); }
  }

void xinit( pthread_cond_t * const cond )
  {
  const int errcode = pthread_cond_init( cond, 0 );
  if( errcode )
    { show_error( "pthread_cond_init", errcode ); cleanup_and_fail(); }
  }


void xdestroy( pthread_mutex_t * const mutex )
  {
  const int errcode = pthread_mutex_destroy( mutex );
  if( errcode )
    { show_error( "pthread_mutex_destroy", errcode ); cleanup_and_fail(); }
  }

void xdestroy( pthread_cond_t * const cond )
  {
  const int errcode = pthread_cond_destroy( cond );
  if( errcode )
    { show_error( "pthread_cond_destroy", errcode ); cleanup_and_fail(); }
  }


void xlock( pthread_mutex_t * const mutex )
  {
  const int errcode = pthread_mutex_lock( mutex );
  if( errcode )
    { show_error( "pthread_mutex_lock", errcode ); cleanup_and_fail(); }
  }


void xunlock( pthread_mutex_t * const mutex )
  {
  const int errcode = pthread_mutex_unlock( mutex );
  if( errcode )
    { show_error( "pthread_mutex_unlock", errcode ); cleanup_and_fail(); }
  }


void xwait( pthread_cond_t * const cond, pthread_mutex_t * const mutex )
  {
  const int errcode = pthread_cond_wait( cond, mutex );
  if( errcode )
    { show_error( "pthread_cond_wait", errcode ); cleanup_and_fail(); }
  }


void xsignal( pthread_cond_t * const cond )
  {
  const int errcode = pthread_cond_signal( cond );
  if( errcode )
    { show_error( "pthread_cond_signal", errcode ); cleanup_and_fail(); }
  }


void xbroadcast( pthread_cond_t * const cond )
  {
  const int errcode = pthread_cond_broadcast( cond );
  if( errcode )
    { show_error( "pthread_cond_broadcast", errcode ); cleanup_and_fail(); }
  }


// namespace {

unsigned long long in_size = 0;
unsigned long long out_size = 0;
// const char * const mem_msg = "Not enough memory. Try a smaller dictionary size";


       // split data from input file into chunks and pass them to
       // courier for packaging and distribution to workers.
extern "C" void * csplitter( void * arg )
  {
  const Splitter_arg & tmp = *(Splitter_arg *)arg;
  Packet_courier & courier = *tmp.courier;
  const Pretty_print & pp = *tmp.pp;
  uint8_t * data = tmp.data_chunk;
  int outfd = tmp.outfd;
  // std::cerr << "csplitter data:" << strlen( (const char*)data) << " data_size=" << tmp.data_size << std::endl;
  // const int infd = tmp.infd;
  const int data_size = tmp.data_size;
  const int block_size = tmp.block_size;

  for (int ate = 0; ate < data_size; ate += block_size) {
  // for( bool first_post = true; ; first_post = false ) {
    // uint8_t * const data = new( std::nothrow ) uint8_t[data_size]; // the main thread would allocate this data
    // if( !data ) { pp( mem_msg ); cleanup_and_fail(); }
    // const int size = readblock( infd, data, data_size );
    // if( size != data_size && errno )
      // { pp(); show_error( "Read error", errno ); cleanup_and_fail(); }

    // if( size > 0 || first_post )	// first packet may be empty
    //   {
    //   in_size += size;
    //   courier.receive_packet( data, size );
    //   if( size < data_size ) break;	// EOF
    //   }
    // else
    //   {
    //   delete[] data;
    //   break;
    //   }
    // }
    int remainder = data_size - ate;
    int size = std::min(block_size, remainder);
    uint8_t * block_data = new( std::nothrow ) uint8_t[ size ];
    memcpy(block_data, &(data[ate]), size );

    // if( size > 0 ) { // first packet may be empty
      // std::cerr << "passing a packet to courier " << std::endl;
      in_size += size;
      courier.receive_packet( block_data, size, outfd );
      // if( size < data_size ) break; // EOF
      // }
    // else {
      // delete[] data;
      // break;
      // }
    }
  // TODO: clear the vector
  // can clear in the caller, but this way memory is freed as soon as all the packets were sent out
    // data->clear();
  // delete [] data;

  courier.finish();			// no more packets to send
  return 0;
  }



       // get packets from courier, replace their contents, and return
       // them to courier.
extern "C" void * cworker( void * arg )
  {
  // std::cerr << "cworker: launched, hanging out" << std::endl;
  const Worker_arg & tmp = *(Worker_arg *)arg;
  Packet_courier & courier = *tmp.courier;
  // const Pretty_print & pp = *tmp.pp;
  const int dictionary_size = tmp.dictionary_size;
  const int match_len_limit = tmp.match_len_limit;

  while( true ) {
    Packet * const packet = courier.distribute_packet();
    if( !packet ) { 
      // std::cerr << "worker::no more packets, finishing " << pthread_self() << std::endl;
      break;		// no more packets to process
    }

    // std::cerr << "got a packet!" << std::endl;

    const int max_compr_size = 42 + packet->size + ( ( packet->size + 7 ) / 8 );
    uint8_t * const new_data = new( std::nothrow ) uint8_t[max_compr_size];
    if( !new_data ) { 
      std::cerr << "FAIL1" << std::endl;
      /*pp( mem_msg );*/ cleanup_and_fail(); 
    }
    const int dict_size = std::max( LZ_min_dictionary_size(),
                                    std::min( dictionary_size, packet->size ) );
    LZ_Encoder * const encoder =
      LZ_compress_open( dict_size, match_len_limit, LLONG_MAX );
    if( !encoder || LZ_compress_errno( encoder ) != LZ_ok )
      {
      std::cerr << "FAIL2" << std::endl;
      if( !encoder || LZ_compress_errno( encoder ) == LZ_mem_error ) {
        // pp( mem_msg );

      }
      else
        internal_error( "invalid argument to encoder" );
      cleanup_and_fail();
      }

    int written = 0;
    int new_pos = 0;
    while( true ) {
      if( LZ_compress_write_size( encoder ) > 0 ) {
        if( written < packet->size ) {
          const int wr = LZ_compress_write( encoder, packet->data + written,
                                          packet->size - written );
          if( wr < 0 ) internal_error( "library error (LZ_compress_write)" );
          written += wr;
        }
        if( written >= packet->size ) { 
          // std::cerr << "WRITTEN=" << written << std::endl;
          delete[] packet->data; 
          LZ_compress_finish( encoder ); 
        }
      }
      const int rd = LZ_compress_read( encoder, new_data + new_pos,
                                 max_compr_size - new_pos );
      if( rd < 0 ) {
        std::cerr << "FAIL3" << std::endl;
        // pp();
        // if( verbosity >= 0 )
        //   std::fprintf( stderr, "LZ_compress_read error: %s.\n",
        //                 LZ_strerror( LZ_compress_errno( encoder ) ) );
        cleanup_and_fail();
      }
      new_pos += rd;
      if( new_pos > max_compr_size )
        internal_error( "packet size exceeded in worker" );
      if( LZ_compress_finished( encoder ) == 1 ) {
        // std::cerr << "worker::finished compressing a block" << std::endl;
        break;
      }
    }

    if( LZ_compress_close( encoder ) < 0 ) { 
      /*pp( "LZ_compress_close failed" );*/ cleanup_and_fail(); 
    }

    // if( verbosity >= 2 && packet->size > 0 ) show_progress( packet->size );
    packet->data = new_data;
    packet->size = new_pos;
    courier.collect_packet( packet );
    }
  return 0;
}


     // get from courier the processed and sorted packets, and write
     // their contents to the output file.
void muxer( Packet_courier & courier /*, const Pretty_print & pp*/)
  {
  std::vector< const Packet * > packet_vector;
  while( true )
    {
    courier.deliver_packets( packet_vector );
    if( packet_vector.size() == 0 ) {
      // std::cerr << "muxer::do not see any more packets. breaking" << std::endl;
      break;		// all workers exited
    }

    for( unsigned i = 0; i < packet_vector.size(); ++i )
      {
      const Packet * const opacket = packet_vector[i];
      out_size += opacket->size;

      int outfd = opacket->outfd;

      if( outfd >= 0 ) {
        // std::cerr << "muxer::write out a compressed block (size=" << opacket->size << ") to fd=" << outfd << std::endl;
        const int wr = writeblock( outfd, opacket->data, opacket->size );
        if( wr != opacket->size )
          { 
            // std::cerr << "wr=" << wr << " packet_size=" << opacket->size << std::endl;
          /*pp();*/ show_error( "Write error", errno ); cleanup_and_fail(); }
      }
      else {
        std::cerr << "ZZZ" << std::endl;
      }
      // std::cerr << "removing temp compressed data" << std::endl;
      delete[] opacket->data;
      delete opacket;
      }
    }
    // std::cerr << "muxer exited" << std::endl;
  }

// } // end namespace


    // init the courier, then start the splitter and the workers and
    // call the muxer.
int plzip_compress( const int data_size, const int dictionary_size,
              const int match_len_limit, const int num_workers,
              uint8_t * data_chunk,
              /*const int infd,*/ const int outfd,
              const Pretty_print & pp, const int debug_level )
  {
  const int slots_per_worker = 2;
  const int num_slots =
    ( ( num_workers > 1 ) ? num_workers * slots_per_worker : 1 );
  in_size = 0;
  out_size = 0;
  Packet_courier courier( num_workers, num_slots );

  Splitter_arg splitter_arg;
  splitter_arg.courier = &courier;
  splitter_arg.pp = &pp;
  splitter_arg.outfd = outfd;
  // splitter_arg.infd = infd;
  splitter_arg.data_chunk = data_chunk;
  splitter_arg.data_size = data_size;
  splitter_arg.block_size = 2 * std::max(data_size, dictionary_size);

  pthread_t splitter_thread;
  int errcode = pthread_create( &splitter_thread, 0, csplitter, &splitter_arg );
  if( errcode )
    { show_error( "Can't create splitter thread", errcode ); cleanup_and_fail(); }

  Worker_arg worker_arg;
  worker_arg.courier = &courier;
  // worker_arg.pp = &pp;
  worker_arg.dictionary_size = dictionary_size;
  worker_arg.match_len_limit = match_len_limit;

  pthread_t * worker_threads = new( std::nothrow ) pthread_t[num_workers];
  if( !worker_threads ) { pp( mem_msg ); cleanup_and_fail(); }
  for( int i = 0; i < num_workers; ++i )
    {
    errcode = pthread_create( worker_threads + i, 0, cworker, &worker_arg );
    if( errcode )
      { show_error( "Can't create worker threads", errcode ); cleanup_and_fail(); }
    }

  muxer( courier );

  for( int i = num_workers - 1; i >= 0; --i )
    {
    errcode = pthread_join( worker_threads[i], 0 );
    if( errcode )
      { show_error( "Can't join worker threads", errcode ); cleanup_and_fail(); }
    }
  delete[] worker_threads;

  errcode = pthread_join( splitter_thread, 0 );
  if( errcode )
    { show_error( "Can't join splitter thread", errcode ); cleanup_and_fail(); }

  // if( verbosity >= 1 )
  //   {
  //   if( in_size == 0 || out_size == 0 )
  //     std::fprintf( stderr, " no data compressed.\n" );
  //   else
  //     std::fprintf( stderr, "%6.3f:1, %6.3f bits/byte, "
  //                           "%5.2f%% saved, %llu in, %llu out.\n",
  //                   (double)in_size / out_size,
  //                   ( 8.0 * out_size ) / in_size,
  //                   100.0 * ( 1.0 - ( (double)out_size / in_size ) ),
  //                   in_size, out_size );
  //   }

  if( debug_level & 1 )
    std::fprintf( stderr,
      "any worker tried to consume from splitter %8u times\n"
      "any worker had to wait                    %8u times\n"
      "muxer tried to consume from workers       %8u times\n"
      "muxer had to wait                         %8u times\n",
      courier.icheck_counter,
      courier.iwait_counter,
      courier.ocheck_counter,
      courier.owait_counter );

  if( !courier.finished() ) internal_error( "courier not finished" );
  return 0;
  }
