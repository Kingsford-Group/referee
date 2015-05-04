#ifndef QUALITY_CLUSTER_HPP
#define QUALITY_CLUSTER_HPP

#include <memory>
#include <fstream>
#include <string>

#include <compress.h>
#include "OutputBuffer.hpp"

using namespace std;

////////////////////////////////////////////////////////////////
//
// cluster construct that holds similar vectors
//
////////////////////////////////////////////////////////////////
class QualityCluster {
public:

	QualityCluster(Packet_courier * c, bool p = false): courier(c), is_pile(p) {}

	QualityCluster(Packet_courier * c, string & profile, int K, bool p = false): courier(c), is_pile(p) {
		this->profile = profile;
		// build profile kmers
		profile_kmers = shared_ptr<unordered_map<string,int>>(new unordered_map<string,int>() );
		for (int i = 0; i < profile.size() - K + 1; i++) {
			auto kmer = profile.substr(i, K);
			if (profile_kmers->find(kmer) == profile_kmers->end()) 
				(*profile_kmers)[kmer] = 0;
			(*profile_kmers)[kmer]++;
		}
	}

	~QualityCluster() {
		// cerr << "Cluster " << cluster_id << ": " << total_vectors << " vectors" << endl;
	}

	shared_ptr<unordered_map<string,int>> getProfileKmers() {return profile_kmers;} 

	int getProfileSize() {return profile.size(); }

	int getClusterID() {return cluster_id;}

	void setClusterID(int id) {cluster_id = id;}

	int size () {return data.size(); }

	// add string's core, pref, suffix
	void add(string & q_v, int id, string & prefix, string & suffix) {
		total_vectors++;
		data.push_back(q_v);
		ids.push_back(id);
		prefices.push_back(prefix);
		suffices.push_back(suffix);
	}

	// add string as is, w/o breaking it into core, pref and suff
	void add(string & q_v, int id) {
		total_vectors++;
		data.push_back(q_v);
		ids.push_back(id);
	}

	////////////////////////////////////////////////////////////////
	void writeCore(string & core) {
		total_vectors++;
		writeString(core, output_str);
	}

	////////////////////////////////////////////////////////////////
	void writePrefix(string & p) {
		writeString(p, prefix_str);
	}

	////////////////////////////////////////////////////////////////
	void writeSuffix(string & s) {
		writeString(s, suffix_str);
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
		for (int j = 0; j < data.size(); j++) {
			writeCore(data[j]);
			if (!is_pile) {
				writePrefix(prefices[j]);
				writeSuffix(suffices[j]);
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
				// TODO: this is O(n), should use deque
				// TODO: do we need to insert the IDs themselves?
				ids.insert(ids.begin() + i, cluster->ids[j]);
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
	void openOutputStream(string & fname, int K_c) {
		string fname_prefix = fname + ".k=" + to_string(K_c);
		// TODO: set dictionary size, match len parameters
		if (cluster_id < 0 || is_pile) {
			output_str = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (fname_prefix + ".other.lz").c_str() ) );
		}
		else {
			output_str = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (fname_prefix + "." + to_string(cluster_id) + ".lz" ).c_str() ) );
			// prefices, suffices
			prefix_str = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (fname_prefix + "." + to_string(cluster_id) + ".prefix.lz").c_str() ) );
			suffix_str = shared_ptr<OutputBuffer>(new OutputBuffer(courier, (fname_prefix + "." + to_string(cluster_id) + ".suffix.lz").c_str() ) );
		}
	}

	////////////////////////////////////////////////////////////////
	void closeOutputStream() {
		output_str->flush();
		if (prefix_str != nullptr) prefix_str->flush();
		if (suffix_str != nullptr) suffix_str->flush();
	}

private:
	Packet_courier * courier;
	bool is_pile = false;
	int cluster_id = -1;
	size_t total_vectors = 0;
	string profile;
	shared_ptr<unordered_map<string, int>> profile_kmers;
	vector<string> data;
	vector<string> prefices;
	vector<string> suffices;
	vector<int> ids;

	shared_ptr<OutputBuffer> output_str;
	shared_ptr<OutputBuffer> prefix_str;
	shared_ptr<OutputBuffer> suffix_str;

};

#endif
