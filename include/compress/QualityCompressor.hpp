/* 
Quality values compressor: finds and refines the clusters based on the first 2mln qual vectors,
maintains OutputBuffers for the clusters
*/

#ifndef QUALITY_COMP_H
#define QUALITY_COMP_H


#include "QualityCluster.hpp"

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
float d2(shared_ptr<unordered_map<string, int>> profile_kmers, float const f1, string & q, int const K) {
	unordered_map<string, pair<int,int>> counts;
	float f2 = q.size() - K + 1;
	for (auto p : *profile_kmers) {
		auto kmer = p.first;
		if (counts.find(p.first) == counts.end()) 
			counts[kmer] = make_pair(0,0);
		counts[kmer].first = p.second;
	}
	for (int i = 0; i < q.size() - K + 1; i++) {
		auto kmer = q.substr(i, K);
		if (counts.find(kmer) == counts.end()) 
			counts[kmer] = make_pair(0,0);
		counts[kmer].second++;
	}
	// compare kmer profiles
	float d2_ = 0;
	for (auto kmer_c : counts) {
		auto c_p = kmer_c.second;
		auto delta = (c_p.first/f1 - c_p.second/f2);
		d2_ += delta * delta;
	}
	return d2_;
}


////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
string trim_ends_s(string & q_v, char mode, string & prefix, string & suffix) {
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
	return q_v.substr(i, j-i);
}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
char weighted_mode(string & qual, int & f) {
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
	}

	///////////////////////////////////////////////////////////
	void handleRead(char * qual_read, int len) {
		string s;
		// TODO: is this a documented fact that need to add '!'
		for (auto i = 0; i < len; i++) s += qual_read[i] + '!';
		if (observed_vectors <= bootstrap_size) {
			assignQualityVector(s, len, observed_vectors, true);
			if (observed_vectors == bootstrap_size)
				refineClusters();
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
	void assignQualityVector(string & q_v, int q_v_len, int id, bool allowNewClusters) {
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
			// try to assign to an existing cluster
			bool found = false;
			for (auto clust : clusters) {
				float d = d2(clust->getProfileKmers(), clust->getProfileSize() - K_c + 1, q_v, K_c);
				if (d < 0.05) {
					found = true;
					clust->add(q_v, id, prefix, suffix);
					cluster_membership.push_back(clust->getClusterID());
					break;
				}
			}

			if (!found) {
				if (allowNewClusters) {
					// create a new cluster w/ this guy in it
					shared_ptr<QualityCluster> cluster(new QualityCluster(courier, q_v, K_c));
					cluster->add(q_v, id, prefix, suffix);
					clusters.push_back(cluster);
					if (clusters.size() % 100 == 0) cerr << clusters.size() << 
						" (" << (id - others->size() ) << "|" << others->size() << ") ";
				}
				else {
					// store into the generic pile
					others->add(q_v, id);
					cluster_membership.push_back(generic_pile_id);
				}
			}
		}
	}

	// write a quality value w/o splitting it into prefix, core, suffix
	void write(string & s, shared_ptr<QualityCluster> c) {
		auto cig = to_cigar(s);
		c->writeCore( cig );
	}

	void write(string & q_v, string & prefix, string & suffix,shared_ptr<QualityCluster> clust) {
		auto cig = to_cigar(q_v);
		clust->writeCore( cig );
		//cig = to_cigar(prefix);
		//clust->writePrefix(cig);
		clust->writePrefix(prefix);
		//cig = to_cigar(suffix);
		//clust->writeSuffix(cig);
		clust->writeSuffix(suffix);
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
			// try to assign to an existing cluster
			bool found = false;
			for (auto clust : clusters) {
				float d = d2(clust->getProfileKmers(), clust->getProfileSize() - K_c + 1.0f, q_v, K_c);
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
		for (auto i = 0; i < clusters.size(); i++) {
			auto clust = clusters[i];
			clust->setClusterID(clust_id);
			float percent = clust->size() / (float)observed_vectors;
			if (percent > percent_abundance) { // one percent
				cerr << "Cluster " << clust_id << ": " << (percent * 100) << "% of vectors seen so far" << endl;
				// record cluster membership for items in this cluster
				clust->fillOutClusterMembership(cluster_membership);
				clust_id++; // only increment for the clusters that we will retain
			}
			else {
				remove.push_back(i);
				// reassign cluster's data to the "others" pile
				others->mergeCluster(clust);
				// for (auto v : clust->data) {
				// 	others->data.push_back(v);
				// }
				// cluster membership is in the sh-pile by default -- no need to update
			}
		}
		// merge vectors in the clusters that were too small into a sh*tpile

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
