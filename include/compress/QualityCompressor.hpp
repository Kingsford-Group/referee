/* 
Quality values compressor: finds and refines the clusters based on the first 2mln qual vectors,
maintains OutputBuffers for the clusters
*/

#ifndef QUALITY_COMP_H
#define QUALITY_COMP_H

#include <chrono>

#include "QualityCluster.hpp"




chrono::duration<double> elapsed_seconds_d2_loop2;
chrono::duration<double> elapsed_seconds_d2_loop3;


float d2_fast(shared_ptr<unordered_map<int, int>> profile_kmers, float const f1, 
	unordered_map<int,int> const & counts, float const f2) {
	chrono::time_point<std::chrono::system_clock> start = chrono::system_clock::now();
	// compare kmer profiles
	float d2_ = 0;
	// all things in cluster profile, but not in this string AND things shared
	for (auto p : *profile_kmers) {
		auto it = counts.find(p.first);
		if (it == counts.end() )
			d2_ += (p.second / f1) * (p.second / f1);
		else {
			float delta = p.second / f1 - (*it).second / f2;
			d2_ += delta * delta;
		}
	}
	chrono::time_point<std::chrono::system_clock> end = chrono::system_clock::now();
	elapsed_seconds_d2_loop2 += (end - start);

	// counts all things in this string, but not in cluster profile
	start = chrono::system_clock::now();
	for (auto p : counts) {
		if (profile_kmers->find(p.first) == profile_kmers->end()) {
			d2_ += p.second / f2 * p.second/f2;
		}
		// else -- already counted above
	}
	end = chrono::system_clock::now();
	elapsed_seconds_d2_loop3 += (end - start);
	return d2_;
}


////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
chrono::duration<double> elapsed_trim_ends;
string trim_ends_s(string & q_v, char mode, string & prefix, string & suffix) {
	chrono::time_point<std::chrono::system_clock> start = chrono::system_clock::now();
	int i = 0, j = q_v.size() - 1;
	while (i < q_v.size() - 1 ) {
		if (q_v[i] != mode && q_v[i+1] != mode) i++; // stricter
		else {
			i++;
			break;
		}
	}
	while (j > (i + 1) ) {
		if (q_v[j] != mode && q_v[j-1] != mode) j--;
		else {
			// j--;
			break;
		}
	}
	if (i < 8) {
		// do not cut prefix
		prefix = "";
		i = 0;
	}
	else {
		prefix = q_v.substr(0, i);
	}
	if (q_v.size() - j < 8) {
		// do not cut suffix
		suffix = "";
		j = string::npos;
	}
	else {
		suffix = q_v.substr(j);
	}
	chrono::time_point<std::chrono::system_clock> end = chrono::system_clock::now();
	elapsed_trim_ends += (end - start);
	return q_v.substr(i, j-i);
}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
chrono::duration<double> elapsed_mode_time;
char weighted_mode(string & qual, int & f) {
	chrono::time_point<std::chrono::system_clock> start = chrono::system_clock::now();
	// TODO: a more careful implementation: need 90 or less? space or smthing else?
	float weights[90];
	int frequencies[90];
	for (int i = 0; i < 90; i++)  {
		weights[i] = 0;
		frequencies[i] = 0;
	}
	for (auto c: qual) {
		int w = c - ' ';
		weights[c - ' '] += 0.1 * w; // w * w
		frequencies[c-' ']++;
	}
	// find max
	float max_f = 0, max_i = -1;
	for (int i = 0; i < 90; i++) {
		if (weights[i] >= max_f) {
			// max_f = weights[i];
			max_f = frequencies[i];
			max_i = i;
		}
	}
	f = max_f;
	chrono::time_point<std::chrono::system_clock> end = chrono::system_clock::now();
	elapsed_mode_time += (end - start);
	return max_i + ' ';
}

string to_cigar(string const & s) {
	string cigar = "";
	int off = 1, i = 1;
	for (; i < s.size(); i++) {
		if (s[i] == s[i-1]) off++;
		else {
			if (off > 1)
				cigar += to_string(off);
			cigar += s[i-1];
			off = 1;
		}
	}
	if (off > 1)
		cigar += to_string(off);
	cigar += s[i-1];
	return cigar;
}


////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
class QualityCompressor {
public:
	///////////////////////////////////////////////////////////
	QualityCompressor(Packet_courier * c, string const & fname, float pa, int bs, int k):
		courier(c),
		fname(fname),
		percent_abundance(pa), 
		bootstrap_size(bs), 
		K_c(k) {
			others = shared_ptr<QualityCluster>(new QualityCluster(courier, true));
	}

	~QualityCompressor() {
		for (auto cluster : clusters) {
			cluster->closeOutputStream();
		}
		others->closeOutputStream();

		// write out cluster membership
		ofstream cluster_out(fname + ".k=" + to_string(K_c) + ".membership" );
		for (auto m : cluster_membership) cluster_out << m << " ";
		cluster_out.close();

		cerr << endl << "Total times in d2 loops: " << elapsed_seconds_d2_loop1.count() << "," <<
			elapsed_seconds_d2_loop2.count() << "," << elapsed_seconds_d2_loop3.count() << "sec" << endl;
		cerr << "Total mode time: " << elapsed_mode_time.count() << "s" << endl;
		cerr << "Total trimming time: " << elapsed_trim_ends.count() << "s" << endl;
		cerr << "Total writing time: " << elapsed_writes.count() << "s" << endl;
	}

	///////////////////////////////////////////////////////////
	void handleRead(char * qual_read, int len) {
		string s;
		// TODO: is this a documented fact that need to add '!'
		for (auto i = 0; i < len; i++) s += qual_read[i] + '!';
		if (observed_vectors <= bootstrap_size) {
			assignQualityVector(s, len, observed_vectors);
			if (observed_vectors == bootstrap_size) {
				cerr << "Time in d2: loop1 " << elapsed_seconds_d2_loop1.count() << "s" << endl;
				cerr << "Time in d2: loop2 " << elapsed_seconds_d2_loop2.count() << "s" << endl;
				cerr << "Time in d2: loop3 " << elapsed_seconds_d2_loop3.count() << "s" << endl;
				// exit(1);
				refineClusters();
			}
		}
		else {
			// have fixed clusters -- assign read to one of the existing clusters, output
			// assignQualityVector(s, len, observed_vectors, false);
			writeToCluster(s, len, observed_vectors);
		}
		observed_vectors++;
	}

	///////////////////////////////////////////////////////////
	void flush() {
		for (auto c : clusters) c->flush();
		others->flush();
	}

private:

	Packet_courier * courier;

	float percent_abundance = 0.3;

	int observed_vectors = 0;

	int bootstrap_size = 0;

	int K_c = 3;

	int generic_pile_id = 0;

	vector<int> cluster_membership;

	// clusters
	vector<shared_ptr<QualityCluster>> clusters;
	// used to hold quality vectors that do not fit anywhere else
	shared_ptr<QualityCluster> others;//(new QualityCluster());

	string fname; // path to the original sam file

	///////////////////////////////////////////////////////////
	// Assumes that q_v is already reverse according to the reverse complement bit
	///////////////////////////////////////////////////////////
	void assignQualityVector(string & q_v, int q_v_len, int id) {
		int mode_frequency = 0;
		char m = weighted_mode(q_v, mode_frequency);
		if (mode_frequency < 0.26 * q_v_len) {
			others->add(q_v, id);
			// record the the vector when into a general pile
			cluster_membership.push_back(generic_pile_id);
			return;
		}
		else {
			string prefix, suffix;
			q_v = trim_ends_s(q_v, m, prefix, suffix);
			unordered_map<int,int> q_kmers = countKmers(q_v, K_c);
			float f2 = q_v.size() - K_c + 1;
			// try to assign to an existing cluster
			bool found = false;
			for (auto clust : clusters) {
				float d = d2_fast(clust->getProfileKmers(), clust->getProfileSize() - K_c + 1.0, q_kmers, f2);
				if (d < 0.05) {
					found = true;
					clust->add(q_v, id, prefix, suffix);
					cluster_membership.push_back(clust->getClusterID());
					break;
				}
			}

			if (!found) {
				// create a new cluster w/ this guy in it
				shared_ptr<QualityCluster> cluster(new QualityCluster(courier, q_v, K_c, m));
				cluster->add(q_v, id, prefix, suffix);
				clusters.push_back(cluster);
				if (clusters.size() % 100 == 0) cerr << clusters.size() << 
					" (" << (id - others->size() ) << "|" << others->size() << ") ";
			}
		}
	}

	// write a quality value w/o splitting it into prefix, core, suffix
	chrono::duration<double> elapsed_writes;
	void write(string & s, shared_ptr<QualityCluster> c) {
		chrono::time_point<std::chrono::system_clock> start = chrono::system_clock::now();
		// auto cig = to_cigar(s);
		c->writeCore( s );
		chrono::time_point<std::chrono::system_clock> end = chrono::system_clock::now();
		elapsed_writes += (end - start);
	}

	void write(string & q_v, string & prefix, string & suffix,shared_ptr<QualityCluster> clust) {
		chrono::time_point<std::chrono::system_clock> start = chrono::system_clock::now();
		// auto cig = to_cigar(q_v);
		clust->writeCore( q_v );
		//cig = to_cigar(prefix);
		//clust->writePrefix(cig);
		clust->writePrefix(prefix);
		//cig = to_cigar(suffix);
		//clust->writeSuffix(cig);
		clust->writeSuffix(suffix);
		chrono::time_point<std::chrono::system_clock> end = chrono::system_clock::now();
		elapsed_writes += (end - start);
	}

	void writeToCluster(string & q_v, int q_v_len, int id) {
		int mode_frequency = 0;
		char m = weighted_mode(q_v, mode_frequency);
		if (mode_frequency < 0.26 * q_v_len) {
			write(q_v, others);
			// record the the vector when into a general pile
			cluster_membership.push_back(generic_pile_id);
			return;
		}
		else {
			string prefix, suffix;
			q_v = trim_ends_s(q_v, m, prefix, suffix);
			auto q_kmers = countKmers(q_v, K_c);
			float f2 = q_v.size() - K_c + 1;
			// try to assign to an existing cluster
			bool found = false;
			for (auto clust : clusters) {
				float d = d2_fast(clust->getProfileKmers(), clust->getProfileSize() - K_c + 1.0f, q_kmers, f2);
				if (d < 0.05) {
					write(q_v, prefix, suffix, clust);
					cluster_membership.push_back(clust->getClusterID());
					return;
				}
			}

			if (!found) {
				write(q_v, others);
				cluster_membership.push_back(generic_pile_id);
			}
		}
	}

	///////////////////////////////////////////////////////////
	//
	///////////////////////////////////////////////////////////
	void refineClusters() {
		// allocate enough for the bootstrap vectors, initialize to 0 to idicate
		// membership in the generic pile
		cluster_membership.resize(observed_vectors, generic_pile_id);
		cerr << "Initial clusters: " << clusters.size() << endl;
		int clust_id = generic_pile_id + 1; // all other cluster IDs will have ids starting with 1
		vector<int> remove;
		cerr << "[INFO] Filtering clusters that are too small..." << endl;
		for (auto i = 0; i < clusters.size(); i++) {
			cerr << i << " ";
			auto clust = clusters[i];
			float percent = clust->size() / (float)observed_vectors;
			cerr << percent << " ";
			if (percent > percent_abundance) { // one percent
				clust->setClusterID(clust_id);
				cerr << "Cluster " << clust_id << ": " << (percent * 100) << "% of vectors seen so far" << endl;
				// record cluster membership for items in this cluster
				clust->fillOutClusterMembership(cluster_membership);
				clust_id++; // only increment for the clusters that we will retain
			}
			else {
				remove.push_back(i);
				// reassign cluster's data to the "others" pile
				cerr << "merging ";
				// TODO: merge
				// others->mergeCluster(clust);
				// for (auto v : clust->data) {
				// 	others->data.push_back(v);
				// }
				// cluster membership is in the sh-pile by default -- no need to update
				cerr << "next ";
			}
		}
		// merge vectors in the clusters that were too small into a sh*tpile
		cerr << "[INFO] merging vectors from rejected clusters..." << endl;

		while (remove.size() > 0) {
			int idx = remove.back();
			remove.pop_back();
			clusters.erase(clusters.begin() + idx);
		}
		cerr << "Clusters remaining: " << (clusters.size() + 1) << endl;
		cerr << "Quality vectors in generic pile: " << others->size() << 
			" (" << (float)others->size() / observed_vectors * 100 << "%)" << endl;

		// write out clusters, empty their containers, keep references to the output buffers around
		for (auto i = 0; i < clusters.size(); i++) {
			clusters[i]->openOutputStream(fname, K_c);
			// write accumualted vectors to the stream
			clusters[i]->flushBootstrapData();
		}

		others->openOutputStream(fname, K_c);
		others->flushBootstrapData();
	}
};

#endif
