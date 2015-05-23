#ifndef OUTPUT_STREAM_LIB_H
#define OUTPUT_STREAM_LIB_H

#include <vector>
#include <queue>

#include <fcntl.h>
#include <unistd.h>

#include <compress.h>

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

	OutputBuffer(Packet_courier * c, const char * fn, int d = 1<<23, int match_len = 36): 
		courier(c),
		dictionary_size(d),
		match_len_limit(match_len) {
		int flags = O_CREAT | O_WRONLY | o_binary;
		// if( force ) flags |= O_TRUNC; else flags |= O_EXCL;
		// offsets stream
		out_fd = open( fn, flags, outfd_mode );
		cerr << "Opened " << fn << " stream for compressed data (fd=" << out_fd << ")" << endl;
	}

	~OutputBuffer() {
		// cerr << "GenOutStream deconstructor, size: " << data.size() << endl;

		// cerr << "GenOutStream finished " << endl;
		// TODO: close here or in the mutex() ?
		cerr << "Closing fd=" << out_fd << ", processed " << total_bytes << " bytes; ";
		close(out_fd);
	}

	int getFD() {return out_fd; }

	////////////////////////////////////////////////////////////////
	// functions for adding offset data
	////////////////////////////////////////////////////////////////
	friend void addReference(int ref_id, bool first, shared_ptr<OutputBuffer> o_str);

	friend void addOffset(int offset, shared_ptr<OutputBuffer> o_str);

	friend void addOffsetPair(int offset, int occur, shared_ptr<OutputBuffer> o_str);

	friend void addUnsignedByte(uint8_t c, shared_ptr<OutputBuffer> o_str);

	friend void writeOp(const edit_pair & edit, int prev_edit_offset, shared_ptr<OutputBuffer> o_str);

	friend void writeOpNoOffset(const edit_pair & edit, shared_ptr<OutputBuffer> o_str);

	friend void writeSpliceOp(const edit_pair & edit, int prev_edit_offset, short splice_len, shared_ptr<OutputBuffer> o_str);

	friend void writeLongSpliceOp(const edit_pair & edit, int prev_edit_offset, int splice_len, shared_ptr<OutputBuffer> o_str);

	friend void writeReadLen(uint8_t read_len, shared_ptr<OutputBuffer> o_str);

	friend void writeQualVector(char * q, int len, shared_ptr<OutputBuffer> o_str);

	friend void writeString(string & s, shared_ptr<OutputBuffer> o_str);

	friend void writeClip(vector<uint8_t> & clip_bytes, int clip_length, shared_ptr<OutputBuffer> o_str);

	friend void writeName(char * read_name, shared_ptr<OutputBuffer> o_str);

	friend void writeFlags(int flags, int mapq, int rnext, int pnext, int tlen, shared_ptr<OutputBuffer> o_str);

	friend void writeOpt(string opt, shared_ptr<OutputBuffer> o_str);

	friend void writeUnaligned(UnalignedRead & read, shared_ptr<OutputBuffer> o_str);

	////////////////////////////////////////////////////////////////
	// stream size
	int size() { return data.size(); }

	void flush() {
		compressAndWriteOut(true);
	}

private:

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
	void compressAndWriteOut(bool flush_all = false) {
		int data_size = data.size();
		total_bytes += data_size;

		// TODO: is this usually one packet?
		// can we avoid memcpy and just give away this data vector to the packet?
		int condition = 0;
		if (flush_all) condition = data_size;
		else condition = data_size - dictionary_size;
		for (int ate = 0; ate < condition; ate += dictionary_size) {
			int remainder = data_size - ate;
			int size = std::min(dictionary_size, remainder);
			uint8_t * block_data = new( std::nothrow ) uint8_t[ size ];
			// memcpy(block_data, &(data[0]), size );
			for (auto i = 0; i < size; i++) block_data[i] = data[i];
			courier->receive_packet( block_data, size, out_fd ); // associate an output stream to the packet
			// erase the elements that we sent to a packet
			data.erase(data.begin(), data.begin() + size);
		}
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
void addReference(int ref_id, bool first, shared_ptr<OutputBuffer> o_str) {
	if (!first) o_str->data.push_back('\n');
	writeInteger(ref_id, o_str->data);
	o_str->data.push_back(' ');

	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void addOffset(int offset, shared_ptr<OutputBuffer> o_str) {
	for (auto c: to_string(offset))
		o_str->data.push_back(c);
	o_str->data.push_back(' ');
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void addOffsetPair(int offset, int occur, shared_ptr<OutputBuffer> o_str) {	
	for (auto c : to_string(offset))
		o_str->data.push_back(c);
	o_str->data.push_back(':');
	for (auto c : to_string(occur))
		o_str->data.push_back(c);
	o_str->data.push_back(' ');

	if (o_str->timeToDump() ) {
    	o_str->compressAndWriteOut();
    }
}

void addUnsignedByte(uint8_t c, shared_ptr<OutputBuffer> o_str) {
	o_str->data.push_back(c);
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void writeOpNoOffset(const edit_pair & edit, shared_ptr<OutputBuffer> o_str) {
	o_str->data.push_back(edit.edit_op);
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void writeOp(const edit_pair & edit, int prev_edit_offset, shared_ptr<OutputBuffer> o_str) {
	o_str->data.push_back(edit.edit_op);
	o_str->data.push_back(edit.edit_pos - prev_edit_offset);
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void writeSpliceOp(const edit_pair & edit, int prev_edit_offset, short splice_len, shared_ptr<OutputBuffer> o_str) {
	o_str->data.push_back(edit.edit_op);
	o_str->data.push_back(edit.edit_pos - prev_edit_offset);
	// stream << (unsigned char) (n >> 8) << (unsigned char) (n & 255);
	o_str->data.push_back( splice_len >> 8 );
	o_str->data.push_back( splice_len & 255 );
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

// uses 3 bytes to encode the splice length
// operation code E (ascii=69) becomes ascii=197
void writeLongSpliceOp(const edit_pair & edit, int prev_edit_offset, int splice_len, shared_ptr<OutputBuffer> o_str) {
	char op = edit.edit_op | ( (int)1 << 7 );
	o_str->data.push_back(op);
	o_str->data.push_back(edit.edit_pos - prev_edit_offset);
	// stream << (unsigned char) (n >> 8) << (unsigned char) (n & 255);
	o_str->data.push_back( splice_len >> 16 );
	o_str->data.push_back( splice_len >> 8 );
	o_str->data.push_back( splice_len & 255 );
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void writeReadLen(uint8_t read_len, shared_ptr<OutputBuffer> o_str) {
	o_str->data.push_back(read_len);
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void writeQualVector(char * q, int len, shared_ptr<OutputBuffer> o_str) {
	for (auto i = 0; i < len; i++)
		o_str->data.push_back( (uint8_t)(q[i] + '!') );
	o_str->data.push_back('\n');
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void writeClip(vector<uint8_t> & clip_bytes, int clip_length, shared_ptr<OutputBuffer> o_str) {
	for (auto i = 0; i < clip_length; i++) o_str->data.push_back(clip_bytes[i]);
	o_str->data.push_back('\n');
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void writeName(char * read_name, shared_ptr<OutputBuffer> o_str) {
	int len = strlen(read_name);
	for (auto i = 0; i < len; i++)
		o_str->data.push_back(read_name[i]);
	o_str->data.push_back('\n');
	if (o_str->timeToDump() ) o_str->compressAndWriteOut();
}

void writeOpt(string opt, shared_ptr<OutputBuffer> o_str) {
	for (auto c : opt) o_str->data.push_back(c);
	o_str->data.push_back('\n');
	if (o_str->timeToDump()) o_str->compressAndWriteOut();
}

void writeFlags(int flags, int mapq, int rnext, int pnext, int tlen, shared_ptr<OutputBuffer> o_str) {
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
	if (o_str->timeToDump()) o_str->compressAndWriteOut();
}

void writeUnaligned(UnalignedRead & read, shared_ptr<OutputBuffer> o_str) {
	o_str->data.push_back('>');
	// write read id
	for (auto c : read.read_name) o_str->data.push_back(c);
	//int len = strlen(read.read_name);
        //for (auto i = 0; i < len; i++)
        //        o_str->data.push_back(read.read_name[i]);
	o_str->data.push_back('\n');
	// write read seq
	for (auto c : read.seq) o_str->data.push_back(c);
		o_str->data.push_back('\n');
	//write read quals
	for (auto c : read.qual) o_str->data.push_back(c);
	o_str->data.push_back('\n');
	// strand
	o_str->data.push_back('+');
	o_str->data.push_back('\n');
	if (o_str->timeToDump()) o_str->compressAndWriteOut();
}

void writeString(string & s, shared_ptr<OutputBuffer> o_str) {
	for (auto c : s) o_str->data.push_back(c);
	o_str->data.push_back('\n');
	if (o_str->timeToDump()) o_str->compressAndWriteOut();
}

#endif
