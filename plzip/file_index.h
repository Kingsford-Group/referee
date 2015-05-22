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

#ifndef INT64_MAX
#define INT64_MAX  0x7FFFFFFFFFFFFFFFLL
#endif

////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////
class Block {
  long long pos_, size_;		// pos + size <= INT64_MAX

public:
  Block( const long long p, const long long s ) : pos_( p ), size_( s ) {}

  long long pos() const { return pos_; }
  long long size() const { return size_; }
  long long end() const { return pos_ + size_; }

  void pos( const long long p ) { pos_ = p; }
  void size( const long long s ) { size_ = s; }

  bool overlaps( const Block & b ) const
    { return ( pos_ < b.end() && b.pos_ < end() ); }
  void shift( Block & b ) { ++size_; ++b.pos_; --b.size_; }
};


////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////
class File_index {
  struct Member {
    Block dblock, mblock;		// data block, member block

    Member( const long long dp, const long long ds,
            const long long mp, const long long ms )
      : dblock( dp, ds ), mblock( mp, ms ) {}
  };

  std::vector< Member > member_vector;
  std::string error_;
  int retval_;

public:
  File_index( const int infd );

  const std::string & error() const { return error_; }
  int retval() const { return retval_; }

  long long data_end() const
    { if( member_vector.size() ) return member_vector.back().dblock.end();
      else return 0; }

  long long file_end() const
    { if( member_vector.size() ) return member_vector.back().mblock.end();
      else return 0; }

  const Block & dblock( const int i ) const
    { return member_vector[i].dblock; }
  const Block & mblock( const int i ) const
    { return member_vector[i].mblock; }
  int members() const { return (int)member_vector.size(); }
  };
