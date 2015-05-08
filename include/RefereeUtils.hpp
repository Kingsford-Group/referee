#ifndef REFEREE_UTILS_H
#define REFEREE_UTILS_H

// taken from Sailfish::SailfishStringUtils.hpp
uint8_t samToTwoBit[] = {0, /*A*/0, /*C*/1, 0, /*G*/2, 0, 0, 0, /*T*/3,
                                           0, 0, 0, 0, 0, 0, 0};

#include <iostream>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <vector>

#include "FastaReader.h"

using namespace std;

typedef int read_id_t;
typedef int transcript_id_t;
typedef int offset_t;
typedef string read_seq_t;

const int L = 15;
char CODE_ORDER[L+1] = "ACDEHGLNRTVWXYZ";
#define position_t int
#define edit_dist_t vector<unsigned short>

////////////////////////////////////////////////////////////////
//
// count number of reads with edits and edits total
//
////////////////////////////////////////////////////////////////
struct Stats {
  // total stats in terms of inidividual occurrences
  int deletions = 0;
  int insertions = 0;
  int splices = 0;
  int clips = 0;
  int hard_clips = 0;
  int mismatches = 0;
  // stats in terms of reads
  int r_deletions = 0;
  int r_insertions = 0;
  int r_mm = 0;
  int r_splices = 0;
  int r_clips = 0;
  int r_hclips = 0;
  int r_rejected = 0;
  int r_unaligned = 0;
};



void printStats(Stats & stats, int const concordant_edits, int const total_edits, int const cutoff) {
    // cerr << endl;
    cerr << "Reads with: dels = " << stats.r_deletions
         << "\tins = " << stats.r_insertions 
         << "\tsplices  = " << stats.r_splices
         << "\tclips = " << stats.r_clips
         << "\tmm = " << stats.r_mm
         << "\trejected = " << stats.r_rejected << ",\tof them unaligned = "
         << stats.r_unaligned << endl;
         

    // cerr << "Total:\tdels = " << stats.deletions
    //      << "\tins = " << stats.insertions 
    //      << "\tsplices  = " << stats.splices
    //      << "\tclips = " << stats.clips
    //      << "\tmismatches = " << stats.mismatches << endl;

    // cerr << endl << "Positions having " << cutoff << "+ edits, 50% or more consistent: " << concordant_edits 
    //     << " (out of " << total_edits << " total bases with edits)" << endl;
}


inline int getFileLength(ifstream & stream) {
  // get file length
  stream.seekg(0, stream.end);  // move to the end
  int length = stream.tellg();   // get the size in bytes?
  stream.seekg(0, stream.beg);  // move back to the beginning
  return length;
}


////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
size_t getFileSize(ifstream & fff) {
    // get file size
    fff.seekg(0, fff.end);
    size_t length = fff.tellg();
    fff.seekg(0);
    fff.close();
    return length;
}



////////////////////////////////////////////////////////////////
//
// reverse complement a string
//
////////////////////////////////////////////////////////////////
inline void reverse_complement(string & s) {
        reverse(s.begin(), s.end());
        // A -> X
        // T -> A
        // X -> T
        // C -> Y
        // G -> C
        // Y -> G
        for (size_t i = 0; i < s.size(); i++)
                if (s[i] == 'A') s[i] = 'X';
                else if (s[i] == 'C') s[i] = 'Y';
        for (size_t i = 0; i < s.size(); i++) {
                if (s[i] == 'T') s[i] = 'A';
                else if (s[i] == 'G') s[i] = 'C';
        }
        for (size_t i = 0; i < s.size(); i++) {
                if (s[i] == 'X') s[i] = 'T';
                else if (s[i] == 'Y') s[i] = 'G';
        }
}

inline vector<uint8_t> reverse_complement(vector<uint8_t> & v) {
	vector<uint8_t> copy(v);
	reverse(copy.begin(), copy.end() );
	for (int i = 0; i < copy.size(); i++)
		if (copy[i] == 'T') copy[i] = 'A';
		else if (copy[i] == 'A') copy[i] = 'T';
		else if (copy[i] == 'C') copy[i] = 'G';
		else if (copy[i] == 'G') copy[i] = 'C';
	return copy;
}

inline string getMinimizer(vector<uint8_t> & v, int len) {
	vector<string> kmers;
	for (int i = 0; i < v.size() - len + 1; i++)
		kmers.push_back( string(v.begin() + i, v.begin() + i + len) );
	sort(kmers.begin(), kmers.end() );
	return kmers[0];
}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
bool cigar_has_edits(string const & s) {
  return  s.find('I') != string::npos ||
      s.find('D') != string::npos ||
      s.find('N') != string::npos ||
      s.find('S') != string::npos ||
      s.find('H') != string::npos;;
}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
pair<char,int> mode(string & Q) {
  char median;
  unordered_map<char,int> freq;
  for (auto q : Q) {
    if (freq.find(q) == freq.end() )
      freq[q] = 0;
    freq[q]++;
  }
  auto m = max_element(freq.begin(), freq.end(),
    [](pair<char,int> const & a, pair<char,int> const & b){
      return a.second < b.second;
  });
  return *m;
}


#endif
