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
#include <iostream>
#include <string>
#include <vector>
#include <pthread.h>
#include <stdint.h>
#include <unistd.h>
#include <lzlib.h>

#include "lzip.h"
#include "file_index.h"


namespace {

enum { max_packet_size = 1 << 20 };


struct Packet	{	// data block
  uint8_t * data;		// data == 0 means end of member
  int size;			// number of bytes in data (if any)
};


class Packet_courier			// moves packets around
  {
public:
  unsigned ocheck_counter;
  unsigned owait_counter;
private:
  int deliver_worker_id;	// worker queue currently delivering packets
  std::vector< std::queue< Packet * > > opacket_queues;
  int num_working;			// number of workers still running
  const int num_workers;		// number of workers
  const int num_slots;			// max output packets in circulation
  int num_free;				// remaining free output slots
  pthread_mutex_t omutex;
  pthread_cond_t oav_or_exit;	// output packet available or all workers exited
  pthread_cond_t slot_av;		// free output slot available

  Packet_courier( const Packet_courier & );	// declared as private
  void operator=( const Packet_courier & );	// declared as private

public:
  Packet_courier( const int workers, const int slots )
    : ocheck_counter( 0 ), owait_counter( 0 ),
      deliver_worker_id( 0 ),
      opacket_queues( workers ), num_working( workers ),
      num_workers( workers ), num_slots( 8 * slots ), num_free( num_slots )
    { xinit( &omutex ); xinit( &oav_or_exit ); xinit( &slot_av ); }

  ~Packet_courier()
    { xdestroy( &slot_av ); xdestroy( &oav_or_exit ); xdestroy( &omutex ); }

  void worker_finished()
    {
    // notify muxer when last worker exits
    xlock( &omutex );
    if( --num_working == 0 ) xsignal( &oav_or_exit );
    xunlock( &omutex );
    }

  // collect a packet from a worker
  void collect_packet( Packet * const opacket, const int worker_id )
    {
    xlock( &omutex );
    if( opacket->data )
      {
      while( worker_id != deliver_worker_id && num_free <= 0 )
        xwait( &slot_av, &omutex );
      --num_free;
      }
    opacket_queues[worker_id].push( opacket );
    if( worker_id == deliver_worker_id ) xsignal( &oav_or_exit );
    xunlock( &omutex );
    }

  // deliver a packet to muxer
  // if packet data == 0, move to next queue and wait again
  Packet * deliver_packet()
    {
    Packet * opacket = 0;
    xlock( &omutex );
    ++ocheck_counter;
    while( true )
      {
      while( opacket_queues[deliver_worker_id].empty() && num_working > 0 )
        {
        ++owait_counter;
        xwait( &oav_or_exit, &omutex );
        }
      if( opacket_queues[deliver_worker_id].empty() ) break;
      opacket = opacket_queues[deliver_worker_id].front();
      opacket_queues[deliver_worker_id].pop();
      if( opacket->data )
        {
        if( ++num_free == 1 ) xsignal( &slot_av );
        break;
        }
      if( ++deliver_worker_id >= num_workers ) deliver_worker_id = 0;
      xbroadcast( &slot_av );		// restart deliver_worker_id thread
      delete opacket; opacket = 0;
      }
    xunlock( &omutex );
    return opacket;
    }

  bool finished()		// all packets delivered to muxer
    {
    if( num_free != num_slots || num_working != 0 ) return false;
    for( int i = 0; i < num_workers; ++i )
      if( !opacket_queues[i].empty() ) return false;
    return true;
    }
  };


struct Worker_arg
  {
  const File_index * file_index;
  Packet_courier * courier;
  const Pretty_print * pp;
  int worker_id;
  int num_workers;
  int infd;
  };


       // read members from file, decompress their contents, and
       // give the produced packets to courier.
extern "C" void * dworker_o( void * arg )
  {
  const Worker_arg & tmp = *(Worker_arg *)arg;
  const File_index & file_index = *tmp.file_index;
  Packet_courier & courier = *tmp.courier;
  const Pretty_print & pp = *tmp.pp;
  const int worker_id = tmp.worker_id;
  const int num_workers = tmp.num_workers;
  const int infd = tmp.infd;
  const int buffer_size = 65536;

  uint8_t * new_data = new( std::nothrow ) uint8_t[max_packet_size];
  uint8_t * const ibuffer = new( std::nothrow ) uint8_t[buffer_size];
  LZ_Decoder * const decoder = LZ_decompress_open();
  if( !new_data || !ibuffer || !decoder ||
      LZ_decompress_errno( decoder ) != LZ_ok )
    { pp( "Not enough memory" ); cleanup_and_fail(); }
  int new_pos = 0;

  for( int i = worker_id; i < file_index.members(); i += num_workers )
    {
    long long member_pos = file_index.mblock( i ).pos();
    long long member_rest = file_index.mblock( i ).size();

    while( member_rest > 0 )
      {
      while( LZ_decompress_write_size( decoder ) > 0 )
        {
        const int size = std::min( LZ_decompress_write_size( decoder ),
                    (int)std::min( (long long)buffer_size, member_rest ) );
        if( size > 0 )
          {
          if( preadblock( infd, ibuffer, size, member_pos ) != size )
            { pp(); show_error( "Read error", errno ); cleanup_and_fail(); }
          member_pos += size;
          member_rest -= size;
          if( LZ_decompress_write( decoder, ibuffer, size ) != size )
            internal_error( "library error (LZ_decompress_write)" );
          }
        if( member_rest <= 0 ) { LZ_decompress_finish( decoder ); break; }
        }
      while( true )			// read and pack decompressed data
        {
        const int rd = LZ_decompress_read( decoder, new_data + new_pos,
                                           max_packet_size - new_pos );
        if( rd < 0 )
          cleanup_and_fail( decompress_read_error( decoder, pp, worker_id ) );
        new_pos += rd;
        if( new_pos > max_packet_size )
          internal_error( "opacket size exceeded in worker" );
        if( new_pos == max_packet_size ||
            LZ_decompress_finished( decoder ) == 1 )
          {
          if( new_pos > 0 )			// make data packet
            {
            Packet * opacket = new Packet;
            opacket->data = new_data;
            opacket->size = new_pos;
            courier.collect_packet( opacket, worker_id );
            new_pos = 0;  // reset new_pos
            new_data = new( std::nothrow ) uint8_t[max_packet_size];  // allocate new_data for the next piece of raw data
            if( !new_data ) { pp( "Not enough memory" ); cleanup_and_fail(); }
            }
          if( LZ_decompress_finished( decoder ) == 1 )
            {
            LZ_decompress_reset( decoder );	// prepare for new member
            Packet * opacket = new Packet;	// end of member token
            opacket->data = 0;
            opacket->size = 0;
            courier.collect_packet( opacket, worker_id );
            break;
            }
          }
        if( rd == 0 ) break;
        }
      }
    }

  delete[] ibuffer; delete[] new_data;
  if( LZ_decompress_member_position( decoder ) != 0 )
    { pp( "Error, some data remains in decoder" ); cleanup_and_fail(); }
  if( LZ_decompress_close( decoder ) < 0 )
    { pp( "LZ_decompress_close failed" ); cleanup_and_fail(); }
  courier.worker_finished();
  return 0;
  }


     // get from courier the processed and sorted packets, and write
     // their contents to the output file.
void muxer( Packet_courier & courier, const Pretty_print & pp, const int outfd )
  {
  while( true )
    {
    Packet * opacket = courier.deliver_packet();
    if( !opacket ) break;	// queue is empty. all workers exited

    if( outfd >= 0 )
      {
      const int wr = writeblock( outfd, opacket->data, opacket->size );
      if( wr != opacket->size )
        { pp(); show_error( "Write error", errno ); cleanup_and_fail(); }
      }
    delete[] opacket->data;
    delete opacket;
    }
  }

} // end namespace


    // init the courier, then start the workers and call the muxer.
int dec_stdout( const int num_workers, const int infd, const int outfd,
                const Pretty_print & pp, const int debug_level,
                const File_index & file_index )
  {
  const int slots_per_worker = 2;
  const int num_slots = ( ( INT_MAX / num_workers >= slots_per_worker ) ?
                          num_workers * slots_per_worker : INT_MAX );

  Packet_courier courier( num_workers, num_slots );

  Worker_arg * worker_args = new( std::nothrow ) Worker_arg[num_workers];
  pthread_t * worker_threads = new( std::nothrow ) pthread_t[num_workers];
  if( !worker_args || !worker_threads )
    { pp( "Not enough memory" ); cleanup_and_fail(); }
  for( int i = 0; i < num_workers; ++i )
    {
    worker_args[i].file_index = &file_index;
    worker_args[i].courier = &courier;
    worker_args[i].pp = &pp;
    worker_args[i].worker_id = i;
    worker_args[i].num_workers = num_workers;
    worker_args[i].infd = infd;
    const int errcode =
      pthread_create( &worker_threads[i], 0, dworker_o, &worker_args[i] );
    if( errcode )
      { show_error( "Can't create worker threads", errcode ); cleanup_and_fail(); }
    }

  muxer( courier, pp, outfd );

  // join
  for( int i = num_workers - 1; i >= 0; --i )
    {
    const int errcode = pthread_join( worker_threads[i], 0 );
    if( errcode )
      { show_error( "Can't join worker threads", errcode ); cleanup_and_fail(); }
    }
  delete[] worker_threads;
  delete[] worker_args;

  // const unsigned long long in_size = file_index.file_end();
  // const unsigned long long out_size = file_index.data_end();
  // if( verbosity >= 2 && out_size > 0 && in_size > 0 )
  //   std::fprintf( stderr, "%6.3f:1, %6.3f bits/byte, %5.2f%% saved.  ",
  //                 (double)out_size / in_size,
  //                 ( 8.0 * in_size ) / out_size,
  //                 100.0 * ( 1.0 - ( (double)in_size / out_size ) ) );
  // if( verbosity >= 3 )
  //   std::fprintf( stderr, "decompressed size %9llu, size %9llu.  ",
  //                 out_size, in_size );

  // if( verbosity >= 1 ) std::fprintf( stderr, "done\n" );

  if( debug_level & 1 )
    std::fprintf( stderr,
      "muxer tried to consume from workers       %8u times\n"
      "muxer had to wait                         %8u times\n",
      courier.ocheck_counter,
      courier.owait_counter );

  if( !courier.finished() ) internal_error( "courier not finished" );
  return 0;
  }
