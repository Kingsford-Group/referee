#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cassert>

#include "RefereeCompress.hpp"
#include "RefereeDecompress.hpp"

using namespace std;

////////////////////////////////////////////////////////////////
struct Params {
    bool decompress = false; // true for decompress, false for compress
    string input_file;  // path to the input file
    int threads = 4;    // max number of threads to use
    bool seq_only = false;
    bool discard_secondary_alignments = false;
    string ref_file;    // path to the reference sequence in *.fa format
    string location;
};

////////////////////////////////////////////////////////////////
void printUsage() {
    cerr << "Referee -- separable compression for sequence alignments" << endl;
    cerr << "To compress:" << endl << 
        "\treferee [options] -r reference.fa alignments.sam" << endl;
    cerr << "To decompress:" << endl << 
        "\treferee -d [options] -r reference.fa alignments.sam" << endl;
    cerr << "Options:" << endl;
    cerr << "\t-t=N                 number of threads" << endl;
    cerr << "\t--seqOnly            encode sequencing data only" << endl;
    cerr << "\t--discardSecondary   discard secondary alignments" << endl;
    cerr << "\tview chrK:L-M        retrieve data from interval [L,M) on chromosome K" << endl;
    cerr << "\t-h, --help           this help" << endl;
}


////////////////////////////////////////////////////////////////
Params parseParameters(int argc, char * argv[]) {
    Params p;
    for (int i = 1; i < argc; i++) {
        if ( strcmp(argv[i], "-d") == 0) {
            p.decompress = true;
        }
        else if ( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printUsage();
            exit(0);
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
                cerr << "[ERROR] Missing arguments for view" << endl;
                exit(1);
            }
            p.location = argv[i];
            p.decompress = true;
        }
        else {
            p.input_file = argv[i];
        }
    }
    if (argc < 3) {
        cerr << "Not enough arguments" << endl;
        printUsage();
        exit(1);
    }
    if (p.input_file.size() == 0) {
        cerr << "[ERROR] Missing required argument: <input_file>" << endl;
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
    // input num threads
    int numParseThreads = p.threads;
    // if parameter not provided -- take up all free threads
        numParseThreads = std::min( (long)numParseThreads, std::min( num_online, max_workers ) );

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
        // decompressFile(p.input_file, p.ref_file, fname_out, numParseThreads);
        decompressFileSequential(p.input_file, p.ref_file, fname_out, p.location);
        cerr << "Restored file written to " << fname_out << endl;
    }
    return 0;
}
