/*  Plzip - Parallel compressor compatible with lzip
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

#ifndef PLZIP_LIB_H
#define PLZIP_LIB_H

#include <pthread.h>

enum {
  min_dictionary_bits = 12,
  min_dictionary_size = 1 << min_dictionary_bits,
  max_dictionary_bits = 29,
  max_dictionary_size = 1 << max_dictionary_bits,
  min_member_size = 36 };

const char * const mem_msg = "Not enough memory. Try a smaller dictionary size";

class Pretty_print
  {
  std::string name_;
  const char * const stdin_name;
  unsigned longest_name;
  mutable bool first_post;

public:

  explicit Pretty_print(  )
    : stdin_name( "(stdin)" ), longest_name( 0 ), first_post( false )
    {

    }

  explicit Pretty_print( const std::vector< std::string > & filenames )
    : stdin_name( "(stdin)" ), longest_name( 0 ), first_post( false )
    {
    const unsigned stdin_name_len = std::strlen( stdin_name );
    for( unsigned i = 0; i < filenames.size(); ++i )
      {
      const std::string & s = filenames[i];
      const unsigned len = ( ( s == "-" ) ? stdin_name_len : s.size() );
      if( len > longest_name ) longest_name = len;
      }
    if( longest_name == 0 ) longest_name = stdin_name_len;
    }

  void set_name( const std::string & filename )
    {
    if( filename.size() && filename != "-" ) name_ = filename;
    else name_ = stdin_name;
    first_post = true;
    }

  void reset() const { if( name_.size() ) first_post = true; }
  const char * name() const { return name_.c_str(); }
  void operator()( const char * const msg = 0 ) const 
    {
    // if( verbosity >= 0 )
      // {
      if( first_post )
        {
        first_post = false;
        std::fprintf( stderr, "  %s: ", name_.c_str() );
        for( unsigned i = 0; i < longest_name - name_.size(); ++i )
          std::fprintf( stderr, " " );
        if( !msg ) std::fflush( stderr );
        }
      if( msg ) std::fprintf( stderr, "%s.\n", msg );
      // }
    }

  };


inline int real_bits( unsigned value )
  {
  int bits = 0;
  while( value > 0 ) { value >>= 1; ++bits; }
  return bits;
  }


const uint8_t magic_string[4] = { 0x4C, 0x5A, 0x49, 0x50 };	// "LZIP"

struct File_header
  {
  uint8_t data[6];			// 0-3 magic bytes
					//   4 version
					//   5 coded_dict_size
  enum { size = 6 };

  void set_magic() { std::memcpy( data, magic_string, 4 ); data[4] = 1; }
  bool verify_magic() const
    { return ( std::memcmp( data, magic_string, 4 ) == 0 ); }

  uint8_t version() const { return data[4]; }
  bool verify_version() const { return ( data[4] == 1 ); }

  unsigned dictionary_size() const
    {
    unsigned sz = ( 1 << ( data[5] & 0x1F ) );
    if( sz > min_dictionary_size )
      sz -= ( sz / 16 ) * ( ( data[5] >> 5 ) & 7 );
    return sz;
    }

  bool dictionary_size( const unsigned sz )
    {
    if( sz >= min_dictionary_size && sz <= max_dictionary_size )
      {
      data[5] = real_bits( sz - 1 );
      if( sz > min_dictionary_size )
        {
        const unsigned base_size = 1 << data[5];
        const unsigned wedge = base_size / 16;
        for( int i = 7; i >= 1; --i )
          if( base_size - ( i * wedge ) >= sz )
            { data[5] |= ( i << 5 ); break; }
        }
      return true;
      }
    return false;
    }
  };


struct File_trailer
  {
  uint8_t data[20];	//  0-3  CRC32 of the uncompressed data
			//  4-11 size of the uncompressed data
			// 12-19 member size including header and trailer

  enum { size = 20 };

  unsigned data_crc() const
    {
    unsigned tmp = 0;
    for( int i = 3; i >= 0; --i ) { tmp <<= 8; tmp += data[i]; }
    return tmp;
    }

  void data_crc( unsigned crc )
    { for( int i = 0; i <= 3; ++i ) { data[i] = (uint8_t)crc; crc >>= 8; } }

  unsigned long long data_size() const
    {
    unsigned long long tmp = 0;
    for( int i = 11; i >= 4; --i ) { tmp <<= 8; tmp += data[i]; }
    return tmp;
    }

  void data_size( unsigned long long sz )
    {
    for( int i = 4; i <= 11; ++i ) { data[i] = (uint8_t)sz; sz >>= 8; }
    }

  unsigned long long member_size() const
    {
    unsigned long long tmp = 0;
    for( int i = 19; i >= 12; --i ) { tmp <<= 8; tmp += data[i]; }
    return tmp;
    }

  void member_size( unsigned long long sz )
    {
    for( int i = 12; i <= 19; ++i ) { data[i] = (uint8_t)sz; sz >>= 8; }
    }
  };


// defined in compress.cc
int readblock( const int fd, uint8_t * const buf, const int size );
int writeblock( const int fd, const uint8_t * const buf, const int size );
void xinit( pthread_mutex_t * const mutex );
void xinit( pthread_cond_t * const cond );
void xdestroy( pthread_mutex_t * const mutex );
void xdestroy( pthread_cond_t * const cond );
void xlock( pthread_mutex_t * const mutex );
void xunlock( pthread_mutex_t * const mutex );
void xwait( pthread_cond_t * const cond, pthread_mutex_t * const mutex );
void xsignal( pthread_cond_t * const cond );
void xbroadcast( pthread_cond_t * const cond );
int plzip_compress( const int data_size, const int dictionary_size,
              const int match_len_limit, int num_workers,
              uint8_t * data_chunk,
              /*const int infd,*/ const int outfd,
              const Pretty_print & pp, const int debug_level );

// defined in file_index.cc
class File_index;

// defined in dec_stdout.cc
int dec_stdout( const int num_workers, const int infd, const int outfd,
                const Pretty_print & pp, const int debug_level,
                const File_index & file_index );

// defined in dec_stream.cc
int dec_stream( const int num_workers, const int infd, const int outfd,
                const Pretty_print & pp, const int debug_level,
                const bool testing );

// defined in decompress.cc
int preadblock( const int fd, uint8_t * const buf, const int size,
                const long long pos );
int pwriteblock( const int fd, const uint8_t * const buf, const int size,
                 const long long pos );
int decompress_read_error( struct LZ_Decoder * const decoder,
                           const Pretty_print & pp, const int worker_id );
int decompress( int num_workers, const int infd, const int outfd,
                const Pretty_print & pp, const int debug_level,
                const bool testing, const bool infd_isreg );

// defined in main.cc
// extern int verbosity;
// terminate the program
void cleanup_and_fail( const int retval = 1 );

void show_error( const char * const msg, const int errcode = 0,
                 const bool help = false );

void internal_error( const char * const msg );

void show_progress( const int packet_size,
                    const Pretty_print * const p = 0,
                    const struct stat * const in_statsp = 0 );



class Slot_tally
  {
  const int num_slots;				// total slots
  int num_free;					// remaining free slots
  pthread_mutex_t mutex;
  pthread_cond_t slot_av;			// free slot available

  Slot_tally( const Slot_tally & );		// declared as private
  void operator=( const Slot_tally & );		// declared as private

public:
  explicit Slot_tally( const int slots )
    : num_slots( slots ), num_free( slots )
    { xinit( &mutex ); xinit( &slot_av ); }

  ~Slot_tally() { xdestroy( &slot_av ); xdestroy( &mutex ); }

  bool all_free() { return ( num_free == num_slots ); }

  void get_slot()				// wait for a free slot
    {
    xlock( &mutex );
    while( num_free <= 0 ) xwait( &slot_av, &mutex );
    --num_free;
    xunlock( &mutex );
    }

  void leave_slot()				// return a slot to the tally
    {
    xlock( &mutex );
    if( ++num_free == 1 ) xsignal( &slot_av );	// num_free was 0
    xunlock( &mutex );
    }

  void leave_slots( const int slots )		// return slots to the tally
    {
    xlock( &mutex );
    num_free += slots;
    if( num_free == slots ) xsignal( &slot_av );	// num_free was 0
    xunlock( &mutex );
    }
  };




#endif