#ifndef PLZIP_COMPRESS_LIB_H
#define PLZIP_COMPRESS_LIB_H

#include <queue>

#include "lzip.h"

struct Packet     // data block with a serial number
  {
  unsigned id;      // serial number assigned as received
  uint8_t * data;
  int size;     // number of bytes in data (if any)
  int outfd;    // output stream to which this packet belongs
  };


class Packet_courier      // moves packets around
  {
public:
  unsigned icheck_counter;
  unsigned iwait_counter;
  unsigned ocheck_counter;
  unsigned owait_counter;
private:
  unsigned receive_id;      // id assigned to next packet received
  unsigned deliver_id;      // id of next packet to be delivered
  Slot_tally slot_tally;    // limits the number of input packets
  std::queue< Packet * > packet_queue;
  std::vector< const Packet * > circular_buffer;
  int num_working;      // number of workers still running
  const int num_slots;      // max packets in circulation
  pthread_mutex_t imutex;
  pthread_cond_t iav_or_eof;  // input packet available or splitter done
  pthread_mutex_t omutex;
  pthread_cond_t oav_or_exit; // output packet available or all workers exited
  bool eof;       // splitter done

  Packet_courier( const Packet_courier & ); // declared as private
  void operator=( const Packet_courier & ); // declared as private

public:
  Packet_courier( const int workers, const int slots )
    : icheck_counter( 0 ), iwait_counter( 0 ),
      ocheck_counter( 0 ), owait_counter( 0 ),
      receive_id( 0 ), deliver_id( 0 ),
      slot_tally( slots ), circular_buffer( slots, (Packet *) 0 ),
      num_working( workers ), num_slots( slots ), eof( false )
    {
    xinit( &imutex ); xinit( &iav_or_eof );
    xinit( &omutex ); xinit( &oav_or_exit );
    }

  ~Packet_courier()
    {
    xdestroy( &oav_or_exit ); xdestroy( &omutex );
    xdestroy( &iav_or_eof ); xdestroy( &imutex );
    }

  // make a packet with data received from splitter
  void receive_packet( uint8_t * const data, const int size, const int outfd )
    {
    Packet * const ipacket = new Packet;
    ipacket->id = receive_id++; // ensures packets are process in order of their arrival
    ipacket->data = data;
    ipacket->size = size;
    ipacket->outfd = outfd;
    slot_tally.get_slot();    // wait for a free slot
    xlock( &imutex );
    packet_queue.push( ipacket );
    xsignal( &iav_or_eof );
    xunlock( &imutex );
    }

  // distribute a packet to a worker
  Packet * distribute_packet()
    {
    Packet * ipacket = 0;
    xlock( &imutex );
    ++icheck_counter;
    while( packet_queue.empty() && !eof )
      {
      ++iwait_counter;
      xwait( &iav_or_eof, &imutex );
      }
    if( !packet_queue.empty() )
      {
      ipacket = packet_queue.front();
      packet_queue.pop();
      }
    xunlock( &imutex );
    if( !ipacket )
      {
      // notify muxer when last worker exits
      xlock( &omutex );
      if( --num_working == 0 ) xsignal( &oav_or_exit );
      xunlock( &omutex );
      }
    return ipacket;
    }

  // collect a packet from a worker (contains compress bytes)
  void collect_packet( const Packet * const opacket )
    {
    const int i = opacket->id%num_slots;
    xlock( &omutex );
    // id collision shouldn't happen
    if( circular_buffer[i] != 0 )
      internal_error( "id collision in collect_packet" );
    // merge packet into circular buffer
    circular_buffer[i] = opacket;
    if( opacket->id == deliver_id ) xsignal( &oav_or_exit );
    xunlock( &omutex );
    }

  // deliver packets to muxer
  void deliver_packets( std::vector< const Packet * > & packet_vector )
    {
    xlock( &omutex );
    ++ocheck_counter;
    int i = deliver_id % num_slots;
    while( circular_buffer[i] == 0 && num_working > 0 )
      {
      ++owait_counter;
      xwait( &oav_or_exit, &omutex );
      }
    packet_vector.clear();
    while( true )
      {
      const Packet * const opacket = circular_buffer[i];
      if( !opacket ) break;
      packet_vector.push_back( opacket );
      circular_buffer[i] = 0;
      ++deliver_id;
      i = deliver_id % num_slots;
      }
    xunlock( &omutex );
    if( packet_vector.size() )    // return slots to the tally
      slot_tally.leave_slots( packet_vector.size() );
    }

  void finish()     // splitter has no more packets to send
    {
    std::cerr << "Courier finished." << std::endl;
    xlock( &imutex );
    eof = true;
    xbroadcast( &iav_or_eof );
    xunlock( &imutex );
    }

  bool finished()   // all packets delivered to muxer
    {
    if( !slot_tally.all_free() || !eof || !packet_queue.empty() ||
        num_working != 0 ) return false;
    for( int i = 0; i < num_slots; ++i )
      if( circular_buffer[i] != 0 ) return false;
    return true;
    }
  };


struct Worker_arg
  {
  Packet_courier * courier;
  // const Pretty_print * pp;
  int dictionary_size;
  int match_len_limit;
  };


extern "C" void * cworker( void * arg );

void muxer( Packet_courier & courier/*, const Pretty_print & pp*/);

struct Splitter_arg
  {
  Packet_courier * courier;
  const Pretty_print * pp;
  // int infd;
  // TODO: shared_ptr?
  uint8_t * data_chunk; // pointer to an array of bytes to compress
  int data_size; // size of a chunk of data given to the compressor to work with
  int block_size; // desired size of block into which to break the data and to compress blocks in parallel
  int outfd;
  };


#endif