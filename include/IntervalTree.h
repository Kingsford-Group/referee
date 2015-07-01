#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

template <class T, typename K = std::size_t>
class Interval {
public:
    K start;
    K stop;
    T value;
    Interval(K s, K e, const T& v)
        : start(s)
        , stop(e)
        , value(v)
    { }
};

template <class T, typename K>
K intervalStart(const Interval<T,K>& i) {
    return i.start;
}

template <class T, typename K>
K intervalStop(const Interval<T,K>& i) {
    return i.stop;
}

template <class T, typename K>
  std::ostream& operator<<(std::ostream& out, Interval<T,K>& i) {
    out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class T, typename K = std::size_t>
class IntervalStartSorter {
public:
    bool operator() (const Interval<T,K>& a, const Interval<T,K>& b) {
        return a.start < b.start;
    }
};


using namespace std;

////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
class GenomicCoordinate {
public:
    int chromosome = 0;
    int offset = 0;

    GenomicCoordinate() {}

    GenomicCoordinate(int c, int o): chromosome(c), offset(o) {}
};

////////////////////////////////////////////////////////////////
//
// [start, end) -- genomic interval describing the span of data that is
// compresssed within a given data block
//
////////////////////////////////////////////////////////////////
class TrueGenomicInterval {
public:
    // start genomic coordinate for this block
    GenomicCoordinate start;
    // end genomic coordinate for this block
    GenomicCoordinate end;
    // number of alignments seen so far before the start of this block
    unsigned long num_alignments;

    bool is_aligned;

    // parse the input interval -- start and end chromosome may differ
    TrueGenomicInterval(const string & s, unsigned long const num) {
        num_alignments = num;
        int split = s.find('-');
        string start_s = s.substr(0, split);
        int split2 = start_s.find(':');
        start.chromosome = stoi( start_s.substr(0, split2) );
        start.offset = stoi( start_s.substr(split2+1) );
        string end_s = s.substr(split+1);
        split2 = end_s.find(':');
        end.chromosome = stoi( end_s.substr(0, split2) );
        auto aligned_split = end_s.find(' ', split2+1);
        end.offset = stoi( end_s.substr(split2+1, aligned_split) );
        is_aligned = (end_s.back() == 't');
    }

    void print() {
        cerr << "Interval " << start.chromosome << ":" << start.offset << " to " <<
            end.chromosome << ":" << end.offset << endl;
    }
};

////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
// size_t?
class GenomicInterval: public Interval<int,int> {
public:
    int chromosome; // id for a chromosome or transcript to which reads aligned

    GenomicInterval(): chromosome(-1), Interval(-1,-1,-1) {}

    GenomicInterval(int c, int s, int e): chromosome(c), Interval(s,e,0) {
    }
};


////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
class RawDataInterval : public GenomicInterval {
private:
    bool decompressed = false;

    bool is_aligned_with_transcript_start = false;

    std::vector<uint8_t> data;

public:
    // offset from the begining of the file to where the block begins
    size_t byte_offset;

    // block size, including the header and trailer
    size_t block_size;

    size_t num_alignments;

    // expected size of the decompressed data -- available in the block trailer
    size_t decompressed_size;

    RawDataInterval(size_t off, size_t size, size_t ex, int chr, int s, int e, unsigned long num, bool is_aligned):
        byte_offset(off), 
        block_size(size),
        decompressed_size(ex),
        GenomicInterval(chr, s, e),
        num_alignments(num),
        is_aligned_with_transcript_start(is_aligned) {}

    void setDecompressedData(std::vector<uint8_t> & d) {
        data = d;
        decompressed = true;
    }

    bool isAlignedWithTranscriptStart() {
        return is_aligned_with_transcript_start;
    }

    std::vector<uint8_t> getData() {
        if (!decompressed) return std::vector<uint8_t>();
        return data;
    }
};



////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////
template <class T, typename K = std::size_t>
class IntervalTree {

public:
    // typedef Interval<T,K> interval;
    typedef RawDataInterval interval;
    typedef std::vector<interval> intervalVector;
    typedef IntervalTree<T,K> intervalTree;
    
    intervalVector intervals;
    intervalTree* left;
    intervalTree* right;
    K center;

    ////////////////////////////////////////////////////////////////
    IntervalTree<T,K>(void)
        : left(NULL)
        , right(NULL)
        , center(0)
    { }

    ////////////////////////////////////////////////////////////////
    IntervalTree<T,K>(const intervalTree& other) 
        : left(NULL)
        , right(NULL)
    {
        center = other.center;
        intervals = other.intervals;
        if (other.left) {
            left = new intervalTree(*other.left);
        }
        if (other.right) {
            right = new intervalTree(*other.right);
        }
    }

    ////////////////////////////////////////////////////////////////
    IntervalTree<T,K>& operator=(const intervalTree& other) {
        center = other.center;
        intervals = other.intervals;
        if (other.left) {
            left = new intervalTree(*other.left);
        } else {
            if (left) delete left;
            left = NULL;
        }
        if (other.right) {
            right = new intervalTree(*other.right);
        } else {
            if (right) delete right;
            right = NULL;
        }
        return *this;
    }

    ////////////////////////////////////////////////////////////////
    IntervalTree<T,K>(
            intervalVector& ivals,
            std::size_t depth = 16,
            std::size_t minbucket = 64,
            K leftextent = 0,
            K rightextent = 0,
            std::size_t maxbucket = 512
            )
        : left(NULL)
        , right(NULL)
    {

        --depth;
        IntervalStartSorter<T,K> intervalStartSorter;
        if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
            std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            intervals = ivals;
        } else {
            if (leftextent == 0 && rightextent == 0) {
                // sort intervals by start
              std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            }

            K leftp = 0;
            K rightp = 0;
            K centerp = 0;
            
            if (leftextent || rightextent) {
                leftp = leftextent;
                rightp = rightextent;
            } else {
                leftp = ivals.front().start;
                std::vector<K> stops;
                stops.resize(ivals.size());
                transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
                rightp = *max_element(stops.begin(), stops.end());
            }

            //centerp = ( leftp + rightp ) / 2;
            centerp = ivals.at(ivals.size() / 2).start;
            center = centerp;

            intervalVector lefts;
            intervalVector rights;

            for (typename intervalVector::iterator i = ivals.begin(); i != ivals.end(); ++i) {
                interval& interval = *i;
                if (interval.stop < center) {
                    lefts.push_back(interval);
                } else if (interval.start > center) {
                    rights.push_back(interval);
                } else {
                    intervals.push_back(interval);
                }
            }

            if (!lefts.empty()) {
                left = new intervalTree(lefts, depth, minbucket, leftp, centerp);
            }
            if (!rights.empty()) {
                right = new intervalTree(rights, depth, minbucket, centerp, rightp);
            }
        }
    }

    interval getFirstInterval() {
        if (intervals.size() > 0)
            return intervals.front();
        else return interval(0,0,0,-1,-1,-1, 0, false);
    }

    void findOverlapping(K start, K stop, intervalVector& overlapping) const {
        if (!intervals.empty() && ! (stop < intervals.front().start)) {
            for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
                const interval& interval = *i;
                if (interval.stop >= start && interval.start <= stop) {
                    overlapping.push_back(interval);
                }
            }
        }

        if (left && start <= center) {
            left->findOverlapping(start, stop, overlapping);
        }

        if (right && stop >= center) {
            right->findOverlapping(start, stop, overlapping);
        }

    }

    void findContained(K start, K stop, intervalVector& contained) const {
        if (!intervals.empty() && ! (stop < intervals.front().start)) {
            for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
                const interval& interval = *i;
                if (interval.start >= start && interval.stop <= stop) {
                    contained.push_back(interval);
                }
            }
        }

        if (left && start <= center) {
            left->findContained(start, stop, contained);
        }

        if (right && stop >= center) {
            right->findContained(start, stop, contained);
        }
    }

    ~IntervalTree(void) {
        // traverse the left and right
        // delete them all the way down
        if (left) {
            delete left;
        }
        if (right) {
            delete right;
        }
    }

};

#endif
