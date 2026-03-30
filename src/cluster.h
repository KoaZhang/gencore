#ifndef CLUSTER_H
#define CLUSTER_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "pair.h"
#include "options.h"
#include "htslib/sam.h"
#include <vector>
#include <map>
#include "stats.h"

using namespace std;

struct ClusterProcessResult{
    ClusterProcessResult() {}

    vector<Pair*> consensusPairs;
    StatsDelta preStatsDelta;
    StatsDelta postStatsDelta;
};

class Cluster {
public:
    Cluster(Options* opt);
    ~Cluster();

    void dump();
    void addPair(Pair* pair);
    void addRead(bam1_t* b);

    bool matches(Pair* p);
    ClusterProcessResult clusterByUMI(int umiDiffThreshold, bool crossContig);


    int getLeftRef(){return mPairs[0]->getLeftRef();}
    int getRightRef(){return mPairs[0]->getRightRef();}
    int getLeftPos(){return mPairs[0]->getLeftPos();}
    int getRightPos(){return mPairs[0]->getRightPos();}
    int getTLEN(){return mPairs[0]->getTLEN();}
    Pair::MapType getMapType(){return mPairs[0]->getMapType();}

    static bool test();

private:
    static int umiDiff(const string& umi1, const string& umi2);
    static bool isDuplex(const string& umi1, const string& umi2);
    int duplexMerge(Pair* p1, Pair* p2);
    int duplexMergeBam(bam1_t* b1, bam1_t* b2);
    
public:
    map<string, Pair*> mPairs;
    Options* mOptions;
};

#endif
