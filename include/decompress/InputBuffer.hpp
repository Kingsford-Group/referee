#ifndef INPUT_STREAM_LIB_H
#define INPUT_STREAM_LIB_H

#include <vector>
#include <queue>

#include <fcntl.h>
#include <unistd.h>

using namespace std;

#define END_OF_STREAM -1
#define SUCCESS 1
#define END_OF_TRANS -2


class InputBuffer {
public:
	////////////////////////////////////////////////////////////////
	InputBuffer(string const & fname, int bs = 1<<22):
		buffer_size(bs),
		f_in(fname.c_str(), ifstream::in)  {
			// f_in = open(fname.c_str(), ifstream::in);
			if (!f_in) {
				cerr << "[INFO] Could not open file: " << fname << endl;
			}
	}

	////////////////////////////////////////////////////////////////
	~InputBuffer() {
		f_in.close();
	}

	bool opened() {return f_in.is_open();}

	////////////////////////////////////////////////////////////////
	bool hasMoreBytes() {
		f_in.peek(); // peek -- will set the eof bits if reach the end of file
		return bytes.size() > 0 || f_in.good(); // either have bytes in the buffer or have not reached eof
	}

	////////////////////////////////////////////////////////////////
	uint8_t getNextByte() {
		if (bytes.size() < 1) {
			readMore();
		}
		uint8_t c = bytes.front();
		bytes.pop_front();
		// current_index++;
		return c;
	}

	////////////////////////////////////////////////////////////////
	vector<uint8_t> getNextNBytes(int n) {
		if (bytes.size() < n) {
			readMore();
		}
		vector<uint8_t> local_bytes;
		// for (int i = 0; i < n && i < bytes.size(); i++) {
		while (n > 0 && bytes.size() > 0) {
			local_bytes.push_back(bytes.front());
			bytes.pop_front();
			// current_index++;
			n--;
		}
		return local_bytes;
	}

	void popNBytes(int n) {
		if (bytes.size() < n) readMore();

		while (n > 0) {
			bytes.pop_front();
			// current_index++;
			n--;
		}
	}

	// int current_index = 0;
private:
	ifstream f_in;
	deque<uint8_t> bytes;

	int buffer_size;


	////////////////////////////////////////////////////////////////
	void readMore() {
		// int s = bytes.size();
		// get more memory
		// bytes.resize(s + buffer_size, 0);
		vector<uint8_t> chunk(buffer_size, 0);
		f_in.read( (char *) &chunk[0], buffer_size);
		int actually_read = f_in.gcount();
		if (actually_read < buffer_size) {
			// cerr << "Requested " << buffer_size << ", got " << actually_read << endl;
		}
		for (auto c : chunk)
			bytes.push_back(c);
	}
};

#endif