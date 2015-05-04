#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#include "RefereeCompress.hpp"
#include "decompress/Decompressor.hpp"

using namespace std;

struct Params {
    string mode = "-c"; // or "-d"
    string input_file;  // path to the input file
    int threads = -1;    // max number of threads to use
    string compression_mode;// = "--seq"; // --seq or --seqUniq
    string ref_file;    // path to the reference sequence in *.fa format
};

Params parseParameters(int argc, char * argv[]) {
    if (argc < 3) {
        cerr << "[ERROR]Not enough arguments" << endl;
        exit(1);
    }
    Params p;
    for (int i = 1; i < argc; i++) {
        if ( strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "-d") == 0) {
            p.mode = argv[i];
        }
        else if (strcmp(argv[i], "-t") == 0) {
            i++;
            p.threads = stoi(argv[i]);
        }
        else if ( strcmp(argv[i], "--seq") == 0 || strcmp(argv[i], "--seqUniq") == 0) {
            p.compression_mode = argv[i];
        }
        else if ( strcmp(argv[i], "-r") == 0) {
            i++;
            p.ref_file = argv[i]; // path to the reference sequence
        }
        else {
            p.input_file = argv[i];
        }
    }
    return p;
}


////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
int main(int argc, char * argv []) {
    Params p = parseParameters(argc, argv);

    // determine the number of unoccipied threads
    const long num_online = std::max( 1L, sysconf( _SC_NPROCESSORS_ONLN ) );
    // determine the number of max available threads on this machine
    long max_workers = sysconf( _SC_THREAD_THREADS_MAX );
    if( max_workers < 1 || max_workers > INT_MAX / (int)sizeof (pthread_t) )
        max_workers = INT_MAX / sizeof (pthread_t);

    if (p.mode.compare("-c") == 0) {
        ////////////////////////////////////////////////
        //
        // compress
        //
        ////////////////////////////////////////////////
        cerr << "Compressing " << p.input_file << endl;
        cerr << "Reference genome: " << p.ref_file << endl;

        int numParseThreads = p.threads;
        bool aligned_seq_only = false;
        bool unique_seq_only = false;
        aligned_seq_only = (p.compression_mode.compare("--seq") == 0);
        unique_seq_only = (p.compression_mode.compare("--seqUniq") == 0);
        if (unique_seq_only) aligned_seq_only = false;

        // if parameter not provided -- take up all free threads
        numParseThreads = std::min( (long)numParseThreads, std::min( num_online, max_workers ) );
        cerr << "Threads: " << numParseThreads << endl;


        // Compressor c(file_name, numParseThreads);
        // c.compress();
        compressFile(p.input_file, p.ref_file, numParseThreads, aligned_seq_only, unique_seq_only);
        cerr << endl << "Compressed streams written to " << p.input_file << ".*" << endl;
    }
    else if (p.mode.compare("-d") == 0) {
        ////////////////////////////////////////////////
        //
        // decompress
        //
        ////////////////////////////////////////////////
        string fname_out = p.input_file + ".recovered";
        if (p.ref_file.size() == 0) {
            cerr << "[ERROR] Missing a reference path." << endl;
            exit(1);
        }
        cerr << "Decompressing " << p.input_file << endl;
        cerr << "Reference: " << p.ref_file << endl;
        cerr << "Saving the recovered data to: " << fname_out << endl;
        int numParseThreads = 2;

        // decompressSAM(file_name, file_name + ".restored",ref_path);
        Decompressor d(p.input_file, fname_out, p.ref_file);
        d.decompress();
        cerr << "Restored file written to " << fname_out << endl;
    }
    return 0;
}
