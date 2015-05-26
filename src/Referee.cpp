#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cassert>

#include "RefereeCompress.hpp"
#include "RefereeDecompress.hpp"
// #include "decompress/Decompressor.hpp"

using namespace std;

struct Params {
    // string mode = "-c"; // or "-d"
    bool decompress = false; // true for decompress, false for compress
    string input_file;  // path to the input file
    int threads = -1;    // max number of threads to use
    // string compression_mode;// = "--seq"; // --seq or --seqUniq
    bool seq_only = false;
    bool discard_secondary_alignments = false;
    string ref_file;    // path to the reference sequence in *.fa format
    string location;
};


////////////////////////////////////////////////////////////////
Params parseParameters(int argc, char * argv[]) {
    if (argc < 3) {
        cerr << "[ERROR] Not enough arguments" << endl;
        exit(1);
    }
    Params p;
    for (int i = 1; i < argc; i++) {
        if ( strcmp(argv[i], "-d") == 0) {
            p.decompress = true;
        }
        else if (strcmp(argv[i], "-t") == 0) {
            i++;
            p.threads = stoi(argv[i]);
        }
        else if ( strcmp(argv[i], "--seqOnly") == 0) {
            p.seq_only = true;
        }
        else if (strcmp(argv[i], "--discardSecondary") == 0) {
            p.discard_secondary_alignments = true;
        }
        else if ( strcmp(argv[i], "-r") == 0) {
            i++;
            // TODO: check that next arg exists
            p.ref_file = argv[i]; // path to the reference sequence
        }
        else if (strcmp(argv[i], "view") == 0) {
            i++;
            if (i >= argc) {
                // TODO: check that arg is formatted correctly to indcate location
                cerr << "[ERROR] Missing argument for view" << endl;
                exit(1);
            }
            p.location = argv[i];
            p.decompress = true;
        }
        else {
            p.input_file = argv[i];
        }
    }
    if (p.input_file.size() == 0) {
        cerr << "[ERROR] Missing required argument: -i <input_file>" << endl;
        exit(1);
    }
    return p;
}


////////////////////////////////////////////////////////////////
//
// MAIN
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
    int numParseThreads = p.threads;
    // if parameter not provided -- take up all free threads
        numParseThreads = std::min( (long)numParseThreads, std::min( num_online, max_workers ) );
        cerr << "Threads: " << numParseThreads << endl;

    if (!p.decompress) {
        ////////////////////////////////////////////////
        //
        // compress
        //
        ////////////////////////////////////////////////
        cerr << "Compressing " << p.input_file << endl;
        cerr << "Reference genome: " << p.ref_file << endl;
        compressFile(p.input_file, p.ref_file, numParseThreads, p.seq_only, p.discard_secondary_alignments);
        cerr << endl << "Compressed streams written to " << p.input_file << ".*" << endl;
    }
    else {
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

        // Decompressor d(p.input_file, fname_out, p.ref_file);
        // d.decompress();
        // TODO: pass a region to decompress & stream out if "view" parameter is present
        decompressFile(p.input_file, p.ref_file, fname_out, numParseThreads);
        cerr << "Restored file written to " << fname_out << endl;
    }
    return 0;
}
