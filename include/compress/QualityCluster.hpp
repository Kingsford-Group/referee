#ifndef QUALITY_CLUSTER_HPP
#define QUALITY_CLUSTER_HPP

#include <memory>
#include <fstream>
#include <string>
#include <chrono>

#include <compress.h>
#include "OutputBuffer.hpp"

using namespace std;


////////////////////////////////////////////////////////////////
// convert kmer from base 89 ('z' - '!') to base 10
////////////////////////////////////////////////////////////////
chrono::duration<double> elapsed_seconds_d2_loop1;
int range = 'z' - '!' + 1;
// only for 4-mers
// 5mers would require int64_t
int r_base[] = {range*range*range, range*range, range, 1};

////////////////////////////////////////////////////////////////
unordered_map<int, int> countKmers(string const & q_v, int const K) {
	chrono::time_point<std::chrono::system_clock> start = chrono::system_clock::now();
	unordered_map<int, int> kmers;
	for (int i = 0; i < q_v.size() - K + 1; i++) {
		int kmer_int = 0;
		for (int j = 0; j < K; j++) {
			kmer_int += (q_v[i + j] - '!') * r_base[j];
		}
		if (kmers.find(kmer_int) == kmers.end() )
			kmers[kmer_int] = 0;
		kmers[kmer_int]++;
	}
	chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
	elapsed_seconds_d2_loop1 += (end - start);
	return kmers;
}


////////////////////////////////////////////////////////////////
//
// cluster construct that holds similar vectors
//
////////////////////////////////////////////////////////////////
class QualityCluster {
///////////////////////////////////////
private:
	Packet_courier * courier;
	char mode = 0;
	bool is_pile = false;
	int cluster_id = -1;
	size_t total_vectors = 0;
	string profile;
	shared_ptr<unordered_map<int, int>> profile_kmers;
	vector<string> data;
	vector<string> prefices;
	vector<string> suffices;
	vector<int> ids;	// a vector of qual IDs -- its index in the original SAM file

	shared_ptr<OutputBuffer> output_str;
	shared_ptr<OutputBuffer> prefix_str;
	shared_ptr<OutputBuffer> suffix_str;

///////////////////////////////////////
public:

	QualityCluster(Packet_courier * c, bool p = false): courier(c), is_pile(p) {}

	QualityCluster(Packet_courier * c, string & profile, int K, char m, bool p = false): 
		courier(c), 
		is_pile(p), 
		mode(m) {
		this->profile = profile;
		// build profile kmers
		// store them within the class instead of recomputing every time
		auto kmers = countKmers(profile, K);
		profile_kmers = make_shared<unordered_map<int,int>>( kmers );
	}

	~QualityCluster() {
		// cerr << "Cluster " << cluster_id << ": " << total_vectors << " vectors" << endl;
	}

	shared_ptr<unordered_map<int,int>> getProfileKmers() {return profile_kmers;} 

	int getProfileSize() {return profile.size(); }

	int getClusterID() {return cluster_id;}

	void setClusterID(int id) {cluster_id = id;}

	int size () {return data.size(); }

	////////////////////////////////////////////////////////////////
	// add string's core, pref, suffix
	////////////////////////////////////////////////////////////////
	void add(string & q_v, int id, string & prefix, string & suffix) {
		total_vectors++;
		data.push_back(q_v);
		ids.push_back(id);
		prefices.push_back(prefix);
		suffices.push_back(suffix);
	}

	////////////////////////////////////////////////////////////////
	// add string as is, w/o breaking it into core, pref and suff
	////////////////////////////////////////////////////////////////
	void add(string & q_v, int id) {
		total_vectors++;
		data.push_back(q_v);
		ids.push_back(id);
	}

	////////////////////////////////////////////////////////////////
	/*
	# for i in xrange(n):
	# 	symbol = line[i]
	# 	if symbol.isdigit():
	# 		# remap to a 130 range
	# 		symbol =  chr(130 + int(symbol) )
	# 	if symbol != 'J':
	# 		delta = i - prev_pos - 1
	# 		if delta > 0:
	# 			edit_str += str(delta)
	# 		edit_str += symbol
	# 		prev_pos = i
	*/
	void writeCore(string const & core, GenomicCoordinate & currentCoord) {
		total_vectors++;
		string deltas = "";
		if (mode >= '!') {
			// cerr << "Delta on quals!" << endl;
			// transform core into a string of deltas to the mode
			int prev_pos = 0;
			for (int i = 0; i < core.size(); i++) {
				uint8_t symbol = core[i];
				if (isdigit(symbol))
					symbol = (uint8_t)(symbol + (uint8_t)130);
				if (symbol != mode) {
					int delta = i - prev_pos - 1;
					if (delta > 0)
						deltas += to_string(delta);
					deltas.push_back(symbol);
					prev_pos = i;
				}
			}
			// cerr << deltas << endl;
		}
		else {
			deltas = core;
		}
		writeString(deltas, output_str, currentCoord, 0);
	}

	////////////////////////////////////////////////////////////////
	void writePrefix(string & p, GenomicCoordinate & currentCoord) {
		writeString(p, prefix_str, currentCoord, 0);
	}

	////////////////////////////////////////////////////////////////
	void writeSuffix(string & s, GenomicCoordinate & currentCoord) {
		writeString(s, suffix_str, currentCoord, 0);
	}

	void flush() {
		output_str->flush();
		if (!is_pile) {
			prefix_str->flush();
			suffix_str->flush();
		}
	}

	////////////////////////////////////////////////////////////////
	void flushBootstrapData() {
		// TODO: keep track of the gneomic coordinates
		GenomicCoordinate gc;
		for (int j = 0; j < data.size(); j++) {
			writeCore(data[j], gc);
			if (!is_pile) {
				writePrefix(prefices[j], gc);
				writeSuffix(suffices[j], gc);
			}
			// else: prefix and suffix are not filled out
		}
		data.clear();
		prefices.clear();
		suffices.clear();
	}

	////////////////////////////////////////////////////////////////
	void fillOutClusterMembership(vector<int> & cluster_membership) {
		for (auto i : ids) {
			cluster_membership[i] = cluster_id;
		}
	}

	////////////////////////////////////////////////////////////////
	// merge the data from cluster to "this", taking into account their order
	////////////////////////////////////////////////////////////////
	void mergeCluster(shared_ptr<QualityCluster> cluster) {
		if (is_pile) {
			int i = 0, j = 0;
			while (j < cluster->size() && i < data.size()) {
				while (ids[i] < cluster->ids[j]) {
					i++;
				}
				// insert before i-th element
				// TODO: this is O(n), should use another data structure
				ids.push_back(cluster->ids[j]); // order does not matter -- it is an index
				data.insert( data.begin() + i, cluster->prefices[j] + cluster->data[j] + cluster->suffices[j] );
				j++;
				i++;
			}
			while (j < cluster->size()) {
				ids.push_back(cluster->ids[j]);
				data.push_back( cluster->prefices[j] + cluster->data[j] + cluster->suffices[j] );
				j++;
			}
		}
		// else: should not happen at this iteration
		// TODO: may merge non-empty clusters later
	}

	////////////////////////////////////////////////////////////////
	void openOutputStream(string & fname_prefix, shared_ptr<ofstream> genomic_coords_out, int K_c) {
		if (cluster_id < 0 || is_pile) {
			output_str = shared_ptr<OutputBuffer>(new OutputBuffer(courier, 
				genomic_coords_out, fname_prefix, ".quals.other.lz", 3 << 20,  12 ) ); // equivalent of -4
		}
		else {
			output_str = shared_ptr<OutputBuffer>(new OutputBuffer(courier, 
				genomic_coords_out, fname_prefix, ".quals." + to_string(cluster_id) + ".lz" ) );
			// prefices, suffices
			prefix_str = shared_ptr<OutputBuffer>(new OutputBuffer(courier, 
				genomic_coords_out, fname_prefix, ".quals." + to_string(cluster_id) + ".prefix.lz" ) );
			suffix_str = shared_ptr<OutputBuffer>(new OutputBuffer(courier, 
				genomic_coords_out, fname_prefix, ".quals." + to_string(cluster_id) + ".suffix.lz" ) );
		}
	}

	////////////////////////////////////////////////////////////////
	void closeOutputStream() {
		output_str->flush();
		if (prefix_str != nullptr) prefix_str->flush();
		if (suffix_str != nullptr) suffix_str->flush();
	}
};

#endif
