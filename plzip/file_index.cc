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

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <stdint.h>
#include <unistd.h>

#include "lzip.h"
#include "file_index.h"


int seek_read( const int fd, uint8_t * const buf, const int size,
               const long long pos ) {
  if( lseek( fd, pos, SEEK_SET ) == pos )
    return readblock( fd, buf, size );
  return 0;
}


const char * format_num( unsigned long long num,
                         unsigned long long limit = -1ULL,
                         const int set_prefix = 0 ) {
  const char * const si_prefix[8] = { "k", "M", "G", "T", "P", "E", "Z", "Y" };
  const char * const binary_prefix[8] = { "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi", "Yi" };
  static bool si = true;
  static char buf[32];

  if( set_prefix ) si = ( set_prefix > 0 );
  const unsigned factor = ( si ? 1000 : 1024 );
  const char * const * prefix = ( si ? si_prefix : binary_prefix );
  const char * p = "";
  bool exact = ( num % factor == 0 );

  for( int i = 0; i < 8 && ( num > limit || ( exact && num >= factor ) ); ++i ) { 
      num /= factor; 
      if( num % factor != 0 ) 
        exact = false; 
      p = prefix[i]; 
  }
  snprintf( buf, sizeof buf, "%llu %s", num, p );
  return buf;
}

////////////////////////////////////////////////////////////////
//
// checks that the file is not corrupted and all information is present to decompress it
//
////////////////////////////////////////////////////////////////
File_index::File_index( const int infd ) : retval_( 0 ) {
  // number of bytes in the compressed file
  const long long isize = lseek( infd, 0, SEEK_END ); 
  // checks for size
  if( isize < 0 ) { 
    error_ = "Input file is not seekable :";
    error_ += std::strerror( errno ); 
    retval_ = 1; 
    return; 
  }
  if( isize > INT64_MAX ) { 
    error_ = "Input file is too long (2^63 bytes or more)";
    retval_ = 2; 
    return; 
  }
  // reading position in the input file
  long long pos = isize;		// always points to a header or EOF
  File_header header;
  File_trailer trailer;

  if( isize < min_member_size )
    { error_ = "Input file is too short"; retval_ = 2; return; }
  if( seek_read( infd, header.data, File_header::size, 0 ) != File_header::size )
    { error_ = "Error reading member header :";
      error_ += std::strerror( errno ); retval_ = 1; return; }
  if( !header.verify_magic() )
    { error_ = "Bad magic number (file not in lzip format)";
      retval_ = 2; return; }
  if( !header.verify_version() )
    { error_ = "Version "; error_ += format_num( header.version() );
      error_ += "member format not supported"; retval_ = 2; return; }

  while( pos >= min_member_size ) {
    std::cerr << "reading at pos=" << pos << " ";
    if( seek_read( infd, trailer.data, File_trailer::size, pos - File_trailer::size ) != File_trailer::size ) { 
      error_ = "Error reading member trailer :";
      error_ += std::strerror( errno ); 
      retval_ = 1; 
      break; 
    }

    const long long member_size = trailer.member_size();
    if ( member_size < min_member_size || member_size > pos ) {
      if( member_vector.size() == 0 )	{	// maybe trailing garbage
        --pos; 
        continue; 
      }
      error_ = "Member size in trailer is corrupt at pos ";
      error_ += format_num( pos - 8 ); retval_ = 2; break;
    }

    if ( seek_read( infd, header.data, File_header::size, pos - member_size ) != File_header::size ) { 
      error_ = "Error reading member header :";
      error_ += std::strerror( errno ); 
      retval_ = 1; 
      break;
    }

    if ( !header.verify_magic() || !header.verify_version() ) {
      if( member_vector.size() == 0 )	{	// maybe trailing garbage
        --pos; 
        continue; 
      }
      error_ = "Bad header at pos ";
      error_ += format_num( pos - member_size ); retval_ = 2; break;
    }

    if( member_vector.size() == 0 && isize - pos > File_header::size &&
          seek_read( infd, header.data, File_header::size, pos ) == File_header::size &&
          header.verify_magic() && header.verify_version() ) {
      error_ = "Last member in input file is truncated or corrupt";
      retval_ = 2; break;
    }
    pos -= member_size;
    member_vector.push_back( Member( 0, trailer.data_size(), pos, member_size ) );
  }

  if( pos != 0 || member_vector.size() == 0 ) {
    member_vector.clear();
    if( retval_ == 0 ) { error_ = "Can't create file index"; retval_ = 2; }
    return;
  }
  std::reverse( member_vector.begin(), member_vector.end() );
  for( unsigned i = 0; i < member_vector.size() - 1; ++i ) {
    const long long end = member_vector[i].dblock.end();
    if( end < 0 || end > INT64_MAX ) {
      member_vector.clear();
      error_ = "Data in input file is too long (2^63 bytes or more)";
      retval_ = 2; return;
    }
    member_vector[i+1].dblock.pos( end );
  }
}
