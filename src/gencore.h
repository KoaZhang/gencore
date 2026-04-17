#ifndef GENCORE_H
#define GENCORE_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "htslib/sam.h"
#include "options.h"
#include "cluster.h"
#include "stats.h"
#include "bed.h"
#include "htslib/sam.h"
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <mutex>
#include <thread>
#include <future>
#include "bamutil.h"

using namespace std;

struct bamComp{
    bool operator()(const bam1_t* b1, const bam1_t* b2) const {
        if(b1->core.tid >= 0) {
            if(b2->core.tid<0 )
                return true;
            else if(b2->core.tid >  b1->core.tid)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos >  b1->core.pos)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid >  b1->core.mtid)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid == b1->core.mtid && b2->core.mpos >  b1->core.mpos)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid == b1->core.mtid && b2->core.mpos == b1->core.mpos) {
                if(b2->core.isize > b1->core.isize)
                    return true;
                else if(b2->core.isize == b1->core.isize && (long)b2->data > (long)b1->data) return true;
                else return false;
            } else
                return false;
        } else {
            if(b2->core.tid<0) {
                return (long)b2->data > (long)b1->data;
            }
            else
                return false;
        }
    }
};

struct ContigResult {
    vector<bam1_t*> reads;
    Stats* postStats;
    Stats* preStats;
};

class Gencore {
public:
    Gencore(Options *opt);
    ~Gencore();

    void consensus();

private:
	void releaseClusters(map<int, map<int, unordered_map<long, Cluster*>>>& clusters);
	void dumpClusters(map<int, map<int, unordered_map<long, Cluster*>>>& clusters);
	void addToCluster(bam1_t* b);
	void addToProperCluster(bam1_t* b);
	void addToUnProperCluster(bam1_t* b);
	void createCluster(map<int, map<int, unordered_map<long, Cluster*>>>& clusters, int tid, int left, long right);
    void outputPair(Pair* p);
    bool outputBam(bam1_t* b);
    void finishConsensus(map<int, map<int, unordered_map<long, Cluster*>>>& clusters);
    void report();
    void outputBam(bam1_t* b, bool isLeft);
    void outputOutSet();
    void writeBam(bam1_t* b);

    ContigResult processContig(int tid, vector<bam1_t*>& reads);
    ContigResult processContigChunk(int tid, vector<bam1_t*>& reads, int startIdx, int endIdx);
    void processCompletedClusters(map<int, map<int, unordered_map<long, Cluster*>>>& clusters, int curPos, ContigResult& result);
    void processAllClusters(map<int, map<int, unordered_map<long, Cluster*>>>& clusters, ContigResult& result);
    void consensusParallel();

private:
    string mInput;
    string mOutput;
    Options *mOptions;
    map<int, map<int, unordered_map<long, Cluster*>>> mProperClusters;
    map<int, map<int, unordered_map<long, Cluster*>>> mUnProperClusters;
    bam_hdr_t *mBamHeader;
    samFile* mOutSam;
    Stats* mPreStats;
    Stats* mPostStats;
    set<bam1_t*, bamComp> mOutSet;
    bool mOutSetCleared;
    int mProcessedTid;
    int mProcessedPos;
    bool mProperClustersFinished;
};

#endif
