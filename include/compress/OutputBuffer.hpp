#ifndef OUTPUT_STREAM_LIB_H
#define OUTPUT_STREAM_LIB_H

#include <vector>
#include <queue>

#include <fcntl.h>
#include <unistd.h>

#include <compress.h>

#include "IntervalTree.h"

const mode_t usr_rw = S_IRUSR | S_IWUSR;
const mode_t all_rw = usr_rw | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH;
mode_t outfd_mode = usr_rw;

#ifdef O_BINARY
const int o_binary = O_BINARY;
#else
const int o_binary = 0;
#endif

using namespace std;

////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
class OutputBuffer {
public:

	///////////////////////////////////////////////////////////
	// provides a pointer to a stream for recording genomic coordiantes corresponding
	// to the compressed blocks
	OutputBuffer(Packet_courier * c, shared_ptr<ofstream> genomic_coord_out, 
		string const & fn, string const & suff, int d = 1<<23, int match_len = 36): 
		courier(c),
		genomic_coordinates_out(genomic_coord_out),
		stream_suffix(suff),
		dictionary_size(d),
		match_len_limit(match_len) {
		int flags = O_CREAT | O_WRONLY | o_binary;
		out_fd = open( (fn + suff).c_str(), flags, outfd_mode );
		// cerr << "Opened " << fn << suff << " stream for compressed data (fd=" << out_fd << ")" << endl;
	}

	///////////////////////////////////////////////////////////
	OutputBuffer(Packet_courier * c, string const & fn, string const & suff, 
			int d = 1<<23, int match_len = 36): 
		courier(c),
		stream_suffix(suff),
		dictionary_size(d),
		match_len_limit(match_len) {
		int flags = O_CREAT | O_WRONLY | o_binary;
		out_fd = open( (fn + suff).c_str(), flags, outfd_mode );
		// cerr << "Opened " << fn << suff << " stream for compressed data (fd=" << out_fd << ")" << endl;
	}

	///////////////////////////////////////////////////////////
	~OutputBuffer() {
		// cerr << "Closing fd=" << out_fd << ", processed " << total_bytes << " bytes; ";
		close(out_fd);
	}

	void setInitialCoordinate(int c, int off) {
		startCoord.chromosome = c;
		startCoord.offset = off;
	}

	void setLastCoordinate(int c, int off, size_t num) {
		// cerr << stream_suffix << " last coord: " << c << ":" << off << endl;
		endCoord.chromosome = c;
		endCoord.offset = off;
		num_alignments = num;
	}

	int getFD() {return out_fd; }

	////////////////////////////////////////////////////////////////
	//
	// functions for adding processed data
	//
	////////////////////////////////////////////////////////////////
	friend void addReference(int ref_id, bool first, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void addOffset(int delta, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void addOffsetPair(int delta, int occur, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void addUnsignedByte(uint8_t c, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeBool(bool b, shared_ptr<OutputBuffer> buf, GenomicCoordinate & coord, size_t num);

	friend void writeOp(const edit_pair & edit, int prev_edit_offset, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeOpNoOffset(const edit_pair & edit, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeSpliceOp(const edit_pair & edit, int prev_edit_offset, short splice_len, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeLongSpliceOp(const edit_pair & edit, int prev_edit_offset, int splice_len, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeReadLen(uint8_t read_len, shared_ptr<OutputBuffer> o_str);

	friend void writeQualVector(char * q, int len, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeString(string & s, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeClip(vector<uint8_t> & clip_bytes, int clip_length, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeName(char * read_name, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeFlags(int flags, int mapq, int rnext, int pnext, int tlen, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeOpt(string opt, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num);

	friend void writeUnaligned(UnalignedRead & read, bool seq_only, shared_ptr<OutputBuffer> o_str);

	////////////////////////////////////////////////////////////////
	// stream size
	int size() { return data.size(); }

	void flush() {
		// TODO: last chromosome, max coordinate
		// cerr << "flushing " << stream_suffix << endl;
		GenomicCoordinate g(std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
		compressAndWriteOut(g, num_alignments, true);
	}

private:

	string stream_suffix;

	size_t num_alignments;

	size_t prev_num_alignments = 0;

	GenomicCoordinate startCoord, endCoord;

	shared_ptr<ofstream> genomic_coordinates_out;

	int64_t total_bytes = 0;

	// vector<uint8_t> data;
	deque<uint8_t> data;
	int out_fd; // output file descriptor
	Packet_courier * courier;

	int dictionary_size = 1<<23;

	int match_len_limit = 36; // equivalent to -6 option

	// int packet_size = 2<<22;

	bool timeToDump() {
		return data.size() >= dictionary_size;
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	void compressAndWriteOut(GenomicCoordinate & currentCoord, size_t num_alignments, bool flush_all = false) {
		// cerr << stream_suffix << " data:" << data.size() << " D: " << dictionary_size << endl;
		int data_size = data.size();
		total_bytes += data_size;

		if (!flush_all)
			endCoord = currentCoord;

		// TODO: is this usually one packet?
		// can we avoid memcpy and just give away this data vector to the packet?
		int max_size = 0;
		if (flush_all) 
			max_size = data_size;
		else 
			max_size = (data_size / dictionary_size) * dictionary_size;

		// cerr << "max packet size: " << max_size << endl;

		// coord output is not set for the stream of unaligned reads
		// and chromosome is not set if flushing data out
		if (currentCoord.chromosome != -1 && genomic_coordinates_out != nullptr)
			(*genomic_coordinates_out) << stream_suffix << " " << prev_num_alignments << " " << 
				startCoord.chromosome << ":" << 
				startCoord.offset << "-" <<
				endCoord.chromosome << ":" << endCoord.offset << endl;

		for (int ate = 0; ate < max_size; ate += dictionary_size) {
			// cerr << "sending the block to PLZIP" << endl;
			int remainder = data_size - ate;
			int size = std::min(dictionary_size, remainder);
			uint8_t * block_data = new( std::nothrow ) uint8_t[ size ];
			int i;
			for (i = 0; i < size; i++) block_data[i] = data[i];
			courier->receive_packet( block_data, size, out_fd ); // associate an output stream to the packet
			// erase the elements that we sent to a packet
			// pop_front is better for this
			i = 0;
			while (i < size) {
				i++;
				data.pop_front();
			}
		}
		startCoord = endCoord;
		prev_num_alignments = num_alignments;
	}
};

////////////////////////////////////////////////////////////////
void writeInteger(int n, deque<uint8_t> & out) {
	for (auto c : to_string(n))
		out.push_back(c);
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void addReference(int ref_id, bool first, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num) {
	if (!first) o_str->data.push_back('\n');
	writeInteger(ref_id, o_str->data);
	o_str->data.push_back(' ');

	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

////////////////////////////////////////////////////////////////
void addOffset(int delta, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num) {
	for (auto c: to_string(delta))
		o_str->data.push_back(c);
	o_str->data.push_back(' ');
	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

////////////////////////////////////////////////////////////////
void addOffsetPair(int delta, int occur, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num) {	
	for (auto c : to_string(delta))
		o_str->data.push_back(c);
	o_str->data.push_back(':');
	for (auto c : to_string(occur))
		o_str->data.push_back(c);
	o_str->data.push_back(' ');

	if (o_str->timeToDump() ) {
    	o_str->compressAndWriteOut(coord, num);
    }
}

////////////////////////////////////////////////////////////////
void addUnsignedByte(uint8_t c, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num) {
	o_str->data.push_back(c);
	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

////////////////////////////////////////////////////////////////
void writeBool(bool b, shared_ptr<OutputBuffer> buf, GenomicCoordinate & coord, size_t num) {
	buf->data.push_back(b);
	if (buf->timeToDump() ) buf->compressAndWriteOut(coord, num);
}

void writeOpNoOffset(const edit_pair & edit, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num) {
	o_str->data.push_back(edit.edit_op);
	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

void writeOp(const edit_pair & edit, int prev_edit_offset, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num) {
	o_str->data.push_back(edit.edit_op);
	o_str->data.push_back(edit.edit_pos - prev_edit_offset);
	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

void writeSpliceOp(const edit_pair & edit, int prev_edit_offset, short splice_len, 
		shared_ptr<OutputBuffer> o_str, 
		GenomicCoordinate & coord, size_t num) {
	o_str->data.push_back(edit.edit_op);
	o_str->data.push_back(edit.edit_pos - prev_edit_offset);
	// stream << (unsigned char) (n >> 8) << (unsigned char) (n & 255);
	o_str->data.push_back( splice_len >> 8 );
	o_str->data.push_back( splice_len & 255 );
	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

// uses 3 bytes to encode the splice length
// operation code E (ascii=69) becomes ascii=197
void writeLongSpliceOp(const edit_pair & edit, int prev_edit_offset, 
		int splice_len, 
		shared_ptr<OutputBuffer> o_str, 
		GenomicCoordinate & coord, size_t num) {
	char op = edit.edit_op | ( (int)1 << 7 );
	o_str->data.push_back(op);
	o_str->data.push_back(edit.edit_pos - prev_edit_offset);
	// stream << (unsigned char) (n >> 8) << (unsigned char) (n & 255);
	o_str->data.push_back( splice_len >> 16 );
	o_str->data.push_back( splice_len >> 8 );
	o_str->data.push_back( splice_len & 255 );
	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

void writeReadLen(uint8_t read_len, shared_ptr<OutputBuffer> o_str) {
	o_str->data.push_back(read_len);
	// GenomicCoordinate g;
	// if (o_str->timeToDump() ) o_str->compressAndWriteOut(g, 0);
}

void writeQualVector(char * q, int len, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num) {
	for (auto i = 0; i < len; i++)
		o_str->data.push_back( (uint8_t)(q[i] + '!') );
	o_str->data.push_back('\n');
	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

void writeClip(vector<uint8_t> & clip_bytes, 
	int clip_length, 
	shared_ptr<OutputBuffer> o_str, 
	GenomicCoordinate & coord, size_t num) {
	for (auto i = 0; i < clip_length; i++) o_str->data.push_back(clip_bytes[i]);
	o_str->data.push_back('\n');
	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

void writeName(char * read_name, 
	shared_ptr<OutputBuffer> o_str, 
	GenomicCoordinate & coord, size_t num) {
	int len = strlen(read_name);
	for (auto i = 0; i < len; i++)
		o_str->data.push_back(read_name[i]);
	o_str->data.push_back('\n');
	if (o_str->timeToDump() ) o_str->compressAndWriteOut(coord, num);
}

void writeOpt(string opt, shared_ptr<OutputBuffer> o_str, 
	GenomicCoordinate & coord, size_t num) {
	for (auto c : opt) o_str->data.push_back(c);
	o_str->data.push_back('\n');
	if (o_str->timeToDump()) o_str->compressAndWriteOut(coord, num);
}

void writeFlags(int flags, int mapq, int rnext, int pnext, int tlen, 
	shared_ptr<OutputBuffer> o_str, 
	GenomicCoordinate & coord, size_t num) {
	// o_str->data.push_back( flags >> 8 );
	// o_str->data.push_back( flags & 255 );
	for (auto c : to_string(flags))
		o_str->data.push_back(c);
	o_str->data.push_back(' ');
	for (auto c : to_string(mapq))
		o_str->data.push_back(c);
	o_str->data.push_back(' ');
	for (auto c : to_string(rnext))
		o_str->data.push_back(c);
	o_str->data.push_back(' ');
	for (auto c : to_string(pnext))
		o_str->data.push_back(c);
	o_str->data.push_back(' ');
	for (auto c : to_string(tlen))
		o_str->data.push_back(c);
	o_str->data.push_back('\n');
	if (o_str->timeToDump()) o_str->compressAndWriteOut(coord, num);
}

void writeUnaligned(UnalignedRead & read, bool seq_only, 
	shared_ptr<OutputBuffer> o_str) {
	o_str->data.push_back('>');
	// write read id
	if (!seq_only) {
		for (auto c : read.read_name) o_str->data.push_back(c);
		//int len = strlen(read.read_name);
	        //for (auto i = 0; i < len; i++)
        	//        o_str->data.push_back(read.read_name[i]);
	}
	o_str->data.push_back('\n');
	// write read seq
	for (auto c : read.seq) o_str->data.push_back(c);
	o_str->data.push_back('\n');
	if (!seq_only) {
		//write read quals
		for (auto c : read.qual) o_str->data.push_back(c);
		o_str->data.push_back('\n');
		// strand
		o_str->data.push_back('+');
		o_str->data.push_back('\n');
	}
	GenomicCoordinate g;	// empty coordinate
	if (o_str->timeToDump()) o_str->compressAndWriteOut(g, 0);
}

void writeString(string & s, shared_ptr<OutputBuffer> o_str, GenomicCoordinate & coord, size_t num) {
	for (auto c : s) o_str->data.push_back(c);
	o_str->data.push_back('\n');
	if (o_str->timeToDump()) o_str->compressAndWriteOut(coord, num);
}

#endif
