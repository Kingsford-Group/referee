#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#include "RefereeCompress.hpp"
#include "Decompressor.hpp"
// #include "DecompressSAM.hpp"

using namespace std;


////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
int main(int argc, char * argv []) {
    if (argc < 3) {
        cerr << "Not enough arguments" << endl;
        exit(1);
    }
    string mode = argv[1];

    // determine the number of unoccipied threads
    const long num_online = std::max( 1L, sysconf( _SC_NPROCESSORS_ONLN ) );
    // determine the number of max available threads on this machine
    long max_workers = sysconf( _SC_THREAD_THREADS_MAX );
    if( max_workers < 1 || max_workers > INT_MAX / (int)sizeof (pthread_t) )
        max_workers = INT_MAX / sizeof (pthread_t);

    if (mode.compare("-c") == 0) {
        ////////////////////////////////////////////////
        //
        // compress
        //
        ////////////////////////////////////////////////
        string file_name = argv[2];
        cerr << "Compressing " << file_name << endl;

        int numParseThreads = -1;
        bool aligned_seq_only = false;
        bool unique_seq_only = false;
        if (argc >= 4)
            numParseThreads = stoi(argv[3]);
        if (argc >= 5) {
            string arg4 = argv[4];
            aligned_seq_only = (arg4.compare("--seq") == 0);
            unique_seq_only = (arg4.compare("--seq2") == 0);
            if (unique_seq_only) aligned_seq_only = false;
        }
        // if parameter not provided -- take up all free threads
        numParseThreads = std::min( (long)numParseThreads, std::min( num_online, max_workers ) );
        cerr << "Threads: " << numParseThreads << endl;

        // Compressor c(file_name, numParseThreads);
        // c.compress();
        compressFile(file_name, numParseThreads, aligned_seq_only, unique_seq_only);
    }
    else if (mode.compare("-d") == 0) {
        ////////////////////////////////////////////////
        //
        // decompress
        //
        ////////////////////////////////////////////////
        string file_name = argv[2];
        string fname_out = file_name + ".recovered";
        if (argc < 4) {
            cerr << "[ERROR] Missing a reference path." << endl;
            exit(1);
        }
        string ref_path = argv[3];
        cerr << "Decompressing " << file_name << endl;
        cerr << "Reference: " << ref_path << endl;
        cerr << "Saving the recovered data to: " << fname_out << endl;
        int numParseThreads = 2;

        // decompressSAM(file_name, file_name + ".restored",ref_path);
        Decompressor d(file_name, fname_out, ref_path);
        d.decompress();
    }

    cerr << "Done" << endl;
    return 0;
}
