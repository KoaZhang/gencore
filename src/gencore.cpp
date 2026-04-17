#include "gencore.h"
#include "bamutil.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "reference.h"
#include <limits.h>
#include <algorithm>

Gencore::Gencore(Options *opt){
    mOptions = opt;
    mBamHeader = NULL;
    mOutSam = NULL;
    mPreStats = new Stats(opt);
    mPreStats->setPostStats(false);
    mPostStats = new Stats(opt);
    mPostStats->setPostStats(true);
    mOutSetCleared = false;
    mProcessedTid = -1;
    mProcessedPos = -1;
    mProperClustersFinished = false;
}

Gencore::~Gencore(){
    outputOutSet();
    releaseClusters(mProperClusters);
    releaseClusters(mUnProperClusters);
    if(mBamHeader != NULL) {
        bam_hdr_destroy(mBamHeader);
        mBamHeader = NULL;
    }
    if(mOutSam != NULL) {
        if (sam_close(mOutSam) < 0) {
            cerr << "ERROR: failed to close " << mOutput << endl;
            exit(-1);
        }
    }
    delete mPreStats;
    delete mPostStats;
}

void Gencore::report() {
    JsonReporter jsonreporter(mOptions);
    jsonreporter.report(mPreStats, mPostStats);
    HtmlReporter htmlreporter(mOptions);
    htmlreporter.report(mPreStats, mPostStats);
}

void Gencore::releaseClusters(map<int, map<int, unordered_map<long, Cluster*>>>& clusters) {
    map<int, map<int, unordered_map<long, Cluster*>>>::iterator iter1;
    map<int, unordered_map<long, Cluster*>>::iterator iter2;
    unordered_map<long, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end(); iter1++) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                delete iter3->second;
            }
        }
    }
}

void Gencore::dumpClusters(map<int, map<int, unordered_map<long, Cluster*>>>& clusters) {
    map<int, map<int, unordered_map<long, Cluster*>>>::iterator iter1;
    map<int, unordered_map<long, Cluster*>>::iterator iter2;
    unordered_map<long, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end(); iter1++) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                 iter3->second->dump();
            }
        }
    }
}

void Gencore::outputOutSet() {
    set<bam1_t*, bamComp>::iterator iter;
    for(iter = mOutSet.begin(); iter!=mOutSet.end(); iter++) {
        writeBam(*iter);
        bam_destroy1(*iter);
    }
    mOutSet.clear();
    mOutSetCleared = true;
}

void Gencore::writeBam(bam1_t* b) {
    static int lastTid = -1;
    static int lastPos = -1;
    static bool warnedUnordered = false;
    if(b->core.tid <lastTid || (b->core.tid == lastTid && b->core.pos <lastPos)) {
        if(b->core.tid >=0 && b->core.pos >= 0) {
            if(!warnedUnordered) {
                cerr << "WARNING: The output will be unordered!" << endl;
                warnedUnordered = true;
            }
        }
    }
    if(sam_write1(mOutSam, mBamHeader, b) <0) {
        error_exit("Writing failed, exiting ...");
    }
    lastTid = b->core.tid;
    lastPos = b->core.pos;

    mPostStats->addRead(b);
}

void Gencore::outputBam(bam1_t* b, bool isLeft) {
    pair<set<bam1_t*, bamComp>::iterator,bool> ret = mOutSet.insert(b);
    if(ret.second == false) {
        cerr << "OOPS, found two completely same reads" << endl;
        BamUtil::dump(b);
        BamUtil::dump(*ret.first);
    }
    set<bam1_t*, bamComp>::iterator insertpos = ret.first;
    insertpos++;
    if(isLeft) {
        set<bam1_t*, bamComp>::iterator iter;
        for(iter = mOutSet.begin(); iter!=insertpos; iter++) {
            if(mProcessedPos == -1 || (*iter)->core.tid>mProcessedTid || ((*iter)->core.tid == mProcessedTid && (*iter)->core.pos >= mProcessedPos)) {
                break;
            }
            writeBam(*iter);
            bam_destroy1(*iter);
        }
        mOutSet.erase(mOutSet.begin(), iter);
    }
}

void Gencore::outputPair(Pair* p) {
    mPostStats->addMolecule(1, p->mLeft && p->mRight);

    if(mOutSam == NULL || mBamHeader == NULL)
        return ;

    if(p->mLeft) {
        outputBam(p->mLeft, true);
        p->mLeft =  NULL;
    }
    if(p->mRight) {
        outputBam(p->mRight, false);
        p->mRight =  NULL;
    }
}

ContigResult Gencore::processContig(int tid, vector<bam1_t*>& reads) {
    ContigResult result;
    result.preStats = new Stats(mOptions);
    result.preStats->setPostStats(false);
    result.postStats = new Stats(mOptions);
    result.postStats->setPostStats(true);
    result.postStats->makeGenomeDepthBuf();
    result.postStats->makeBedStats(mPreStats->mBedStats);

    map<int, map<int, unordered_map<long, Cluster*>>> clusters;

    for(int i=0; i<reads.size(); i++) {
        bam1_t* b = reads[i];
        int left = b->core.pos;
        long right;

        if(b->core.mtid == b->core.tid && abs(b->core.mpos - b->core.pos) < 100000) {
            if(b->core.isize < 0) {
                left = b->core.mpos;
            }
            right = left + abs(b->core.isize) - 1;
        } else {
            if(b->core.mtid < 0) {
                result.reads.push_back(b);
                reads[i] = NULL;
                continue;
            } else {
                right = -1L * (long)mBamHeader->target_len[b->core.tid] * (long)(b->core.mtid+1) + (long)b->core.mpos;
            }
        }

        createCluster(clusters, tid, left, right);
        clusters[tid][left][right]->addRead(b);
    }

    map<int, map<int, unordered_map<long, Cluster*>>>::iterator iter1;
    map<int, unordered_map<long, Cluster*>>::iterator iter2;
    unordered_map<long, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end();) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ) {
                if(iter1->first < 0 || iter2->first < 0 ) {
                    unordered_map<string, Pair*>::iterator iterOfPairs;
                    for(iterOfPairs = iter3->second->mPairs.begin(); iterOfPairs!=iter3->second->mPairs.end(); iterOfPairs++) {
                        Pair* p = iterOfPairs->second;
                        result.preStats->addMolecule(1, p->mLeft && p->mRight);
                        result.postStats->addMolecule(1, p->mLeft && p->mRight);
                        if(p->mLeft) {
                            result.reads.push_back(p->mLeft);
                            p->mLeft = NULL;
                        }
                        if(p->mRight) {
                            result.reads.push_back(p->mRight);
                            p->mRight = NULL;
                        }
                        delete p;
                    }
                } else {
                    vector<Pair*> csPairs = iter3->second->clusterByUMI(mOptions->properReadsUmiDiffThreshold, result.preStats, result.postStats, iter3->first < 0);
                    for(int i=0; i<csPairs.size(); i++) {
                        Pair* p = csPairs[i];
                        result.preStats->addMolecule(1, p->mLeft && p->mRight);
                        result.postStats->addMolecule(1, p->mLeft && p->mRight);
                        if(p->mLeft) {
                            result.reads.push_back(p->mLeft);
                            p->mLeft = NULL;
                        }
                        if(p->mRight) {
                            result.reads.push_back(p->mRight);
                            p->mRight = NULL;
                        }
                        delete p;
                    }
                }
                delete iter3->second;
                iter3 = iter2->second.erase(iter3);
            }
            if(iter2->second.size() == 0) {
                iter2 = iter1->second.erase(iter2);
            } else {
                iter2++;
            }
        }
        if(iter1->second.size() == 0) {
            iter1 = clusters.erase(iter1);
        } else {
            iter1++;
        }
    }

    sort(result.reads.begin(), result.reads.end(), bamComp());

    return result;
}

void Gencore::consensusParallel(){
    samFile *in;
    in = sam_open(mOptions->input.c_str(), "r");
    if (!in) {
        cerr << "ERROR: failed to open " << mOptions->input << endl;
        exit(-1);
    }
    hts_set_threads(in, 1);

    if(ends_with(mOptions->output, "sam"))
        mOutSam = sam_open(mOptions->output.c_str(), "w");
    else
        mOutSam = sam_open(mOptions->output.c_str(), "wb");
    if (!mOutSam) {
        cerr << "ERROR: failed to open output " << mOptions->output << endl;
        exit(-1);
    }

    mBamHeader = sam_hdr_read(in);
    mOptions->bamHeader = mBamHeader;
    mPreStats->makeGenomeDepthBuf();
    mPreStats->makeBedStats();
    mPostStats->makeGenomeDepthBuf();
    mPostStats->makeBedStats(mPreStats->mBedStats);

    if (mBamHeader == NULL || mBamHeader->n_targets == 0) {
        cerr << "ERROR: this SAM file has no header " << mInput << endl;
        exit(-1);
    }
    BamUtil::dumpHeader(mBamHeader);

    if (sam_hdr_write(mOutSam, mBamHeader) < 0) {
        cerr << "failed to write header" << endl;
        exit(-1);
    }

    map<int, vector<bam1_t*>> contigReads;
    vector<int> contigOrder;

    bam1_t *b = NULL;
    b = bam_init1();
    int r;
    int count = 0;
    int lastTid = -1;
    int lastPos = -1;
    bool hasPE = false;
    bool isFirst = true;
    while ((r = sam_read1(in, mBamHeader, b)) >= 0) {
        if(isFirst) {
            if(mOptions->umiPrefix == "auto") {
                string umi = BamUtil::getQName(b);
                if(umi.find("umi_") != string::npos)
                    mOptions->umiPrefix = "umi";
                else if(umi.find("UMI_") != string::npos)
                    mOptions->umiPrefix = "UMI";
                else
                    mOptions->umiPrefix = "";

                if(!mOptions->umiPrefix.empty())
                    cerr << endl << "Detected UMI prefix: " << mOptions->umiPrefix << endl << endl;
            }
            isFirst = false;
        }
        mPreStats->addRead(b);
        count++;
        if(count < 1000) {
            if(b->core.mtid >= 0)
                hasPE = true;
        }
        if(count == 1000 && hasPE == false) {
            cerr << "WARNING: seems that the input data is single-end, gencore will not make consensus read and remove duplication for SE data since grouping by coordination will be inaccurate." << endl << endl;
        }

        if(b->core.tid <lastTid || (b->core.tid == lastTid && b->core.pos <lastPos)) {
            if(b->core.tid >=0 && b->core.pos >= 0) {
                cerr << "ERROR: the input is unsorted. Found " << b->core.tid << ":" << b->core.pos << " after " << lastTid << ":" << lastPos << endl;
                cerr << "Please sort the input first." << endl << endl;
                BamUtil::dump(b);
                exit(-1);
            }
        }
        if(mOptions->maxContig>0 && b->core.tid>=mOptions->maxContig){
            b = bam_init1();
            break;
        }
        if(mOptions->debug && b->core.tid > lastTid) {
            cerr << "Starting contig " << b->core.tid << endl;
        }
        lastTid = b->core.tid;
        lastPos = b->core.pos;

        if(b->core.tid < 0 || b->core.pos < 0 ) {
            continue;
        }

        if(!BamUtil::isPrimary(b)) {
            continue;
        }

        int tid = b->core.tid;
        if(contigReads.find(tid) == contigReads.end()) {
            contigOrder.push_back(tid);
        }
        contigReads[tid].push_back(b);
        b = bam_init1();
    }

    bam_destroy1(b);
    sam_close(in);

    int numContigs = contigOrder.size();
    int numThreads = min(mOptions->threads, numContigs);
    cerr << "Processing " << numContigs << " contigs with " << numThreads << " threads" << endl;

    map<int, ContigResult> contigResults;

    if(numThreads <= 1) {
        for(int i=0; i<numContigs; i++) {
            int tid = contigOrder[i];
            if(mOptions->debug) {
                cerr << "Processing contig " << tid << " with " << contigReads[tid].size() << " reads" << endl;
            }
            contigResults[tid] = processContig(tid, contigReads[tid]);
        }
    } else {
        int batchSize = max(1, numContigs / numThreads);
        vector<future<map<int, ContigResult>>> futures;

        for(int t=0; t<numThreads; t++) {
            int startIdx = t * batchSize;
            int endIdx = (t == numThreads - 1) ? numContigs : (t + 1) * batchSize;
            if(startIdx >= numContigs) break;

            vector<int> batchContigs;
            for(int i=startIdx; i<endIdx; i++) {
                batchContigs.push_back(contigOrder[i]);
            }

            futures.push_back(async(launch::async, [this, batchContigs, &contigReads]() {
                map<int, ContigResult> localResults;
                for(int tid : batchContigs) {
                    localResults[tid] = this->processContig(tid, contigReads[tid]);
                }
                return localResults;
            }));
        }

        for(auto& f : futures) {
            map<int, ContigResult> localResults = f.get();
            for(auto& kv : localResults) {
                contigResults[kv.first] = move(kv.second);
            }
        }
    }

    for(int i=0; i<numContigs; i++) {
        int tid = contigOrder[i];
        ContigResult& cr = contigResults[tid];
        for(int j=0; j<cr.reads.size(); j++) {
            if(sam_write1(mOutSam, mBamHeader, cr.reads[j]) < 0) {
                error_exit("Writing failed, exiting ...");
            }
            mPostStats->addRead(cr.reads[j]);
            bam_destroy1(cr.reads[j]);
        }
        mPreStats->merge(cr.preStats);
        mPostStats->merge(cr.postStats);
        delete cr.preStats;
        delete cr.postStats;
    }

    cerr << "----Before gencore processing:" << endl;
    mPreStats->print();

    cerr << endl << "----After gencore processing:" << endl;
    mPostStats->print();

    report();
}

void Gencore::consensusSerial(){
    samFile *in;
    in = sam_open(mOptions->input.c_str(), "r");
    if (!in) {
        cerr << "ERROR: failed to open " << mOptions->input << endl;
        exit(-1);
    }

    if(ends_with(mOptions->output, "sam"))
        mOutSam = sam_open(mOptions->output.c_str(), "w");
    else
        mOutSam = sam_open(mOptions->output.c_str(), "wb");
    if (!mOutSam) {
        cerr << "ERROR: failed to open output " << mOptions->output << endl;
        exit(-1);
    }

    mBamHeader = sam_hdr_read(in);
    mOptions->bamHeader = mBamHeader;
    mPreStats->makeGenomeDepthBuf();
    mPreStats->makeBedStats();
    mPostStats->makeGenomeDepthBuf();
    mPostStats->makeBedStats(mPreStats->mBedStats);

    if (mBamHeader == NULL || mBamHeader->n_targets == 0) {
        cerr << "ERROR: this SAM file has no header " << mInput << endl;
        exit(-1);
    }
    BamUtil::dumpHeader(mBamHeader);

    if (sam_hdr_write(mOutSam, mBamHeader) < 0) {
        cerr << "failed to write header" << endl;
        exit(-1);
    }

    bam1_t *b = NULL;
    b = bam_init1();
    int r;
    int count = 0;
    int lastTid = -1;
    int lastPos = -1;
    bool hasPE = false;
    bool isFirst = true;
    while ((r = sam_read1(in, mBamHeader, b)) >= 0) {
        if(isFirst) {
            if(mOptions->umiPrefix == "auto") {
                string umi = BamUtil::getQName(b);
                if(umi.find("umi_") != string::npos)
                    mOptions->umiPrefix = "umi";
                else if(umi.find("UMI_") != string::npos)
                    mOptions->umiPrefix = "UMI";
                else
                    mOptions->umiPrefix = "";

                if(!mOptions->umiPrefix.empty())
                    cerr << endl << "Detected UMI prefix: " << mOptions->umiPrefix << endl << endl;
            }
            isFirst = false;
        }
        mPreStats->addRead(b);
        count++;
        if(count < 1000) {
            if(b->core.mtid >= 0)
                hasPE = true;
        }
        if(count == 1000 && hasPE == false) {
            cerr << "WARNING: seems that the input data is single-end, gencore will not make consensus read and remove duplication for SE data since grouping by coordination will be inaccurate." << endl << endl;
        }

        if(b->core.tid <lastTid || (b->core.tid == lastTid && b->core.pos <lastPos)) {
            if(b->core.tid >=0 && b->core.pos >= 0) {
                cerr << "ERROR: the input is unsorted. Found " << b->core.tid << ":" << b->core.pos << " after " << lastTid << ":" << lastPos << endl;
                cerr << "Please sort the input first." << endl << endl;
                BamUtil::dump(b);
                exit(-1);
            }
        }
        if(mOptions->maxContig>0 && b->core.tid>=mOptions->maxContig){
            b = bam_init1();
            break;
        }
        if(mOptions->debug && b->core.tid > lastTid) {
            cerr << "Starting contig " << b->core.tid << endl;
        }
        lastTid = b->core.tid;
        lastPos = b->core.pos;

        if(b->core.tid < 0 || b->core.pos < 0 ) {
            if(!mOutSetCleared) {
                if(!mProperClustersFinished) {
                    mProperClustersFinished = true;
                    finishConsensus(mProperClusters);
                }
                outputOutSet();
            }
            continue;
        }

        if(!BamUtil::isPrimary(b)) {
            continue;
        }
        addToCluster(b);
        b = bam_init1();
    }

    if(!mProperClustersFinished) {
        mProperClustersFinished = true;
        finishConsensus(mProperClusters);
    }

    bam_destroy1(b);
    sam_close(in);

    cerr << "----Before gencore processing:" << endl;
    mPreStats->print();

    cerr << endl << "----After gencore processing:" << endl;
    mPostStats->print();

    report();
}

void Gencore::consensus(){
    if(mOptions->threads > 1) {
        consensusParallel();
    } else {
        consensusSerial();
    }
}

void Gencore::addToProperCluster(bam1_t* b) {
    int tid = b->core.tid;
    int left = b->core.pos;
    long right;

    if(b->core.mtid == b->core.tid && abs(b->core.mpos - b->core.pos) < 100000) {
        if(b->core.isize < 0) {
            left = b->core.mpos;
        }
        right = left + abs(b->core.isize) -  1;
    } else {
        if(b->core.mtid < 0) {
            outputBam(b, true);
            return;
        } else {
            right = -1L * (long)mBamHeader->target_len[b->core.tid] * (long)(b->core.mtid+1) + (long)b->core.mpos;
        }
    }

    createCluster(mProperClusters, tid, left, right);
    mProperClusters[tid][left][right]->addRead(b);


    static int tick = 0;
    tick++;
    if(tick % 10000 != 0)
        return;

    map<int, map<int, unordered_map<long, Cluster*>>>::iterator iter1;
    map<int, unordered_map<long, Cluster*>>::iterator iter2;
    unordered_map<long, Cluster*>::iterator iter3;
    bool needBreak = false;
    int curProcessedTid = INT_MAX;
    int curProcessedPos = -1;
    int processedPos;
    for(iter1 = mProperClusters.begin(); iter1 != mProperClusters.end();) {
        if(iter1->first > tid || needBreak) {
            if(curProcessedTid > iter1->first) {
                curProcessedTid = iter1->first;
                curProcessedPos = processedPos;
            }
            break;
        }
        processedPos =mBamHeader->target_len[iter1->first];
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ) {
            if(iter1->first == tid && iter2->first >= b->core.pos) {
                if(processedPos > iter2->first)
                    processedPos = iter2->first;
                needBreak = true;
                break;
            }
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ) {
                if(iter1->first == tid && iter3->first >= b->core.pos) {
                    break;
                }
                vector<Pair*> csPairs = iter3->second->clusterByUMI(mOptions->properReadsUmiDiffThreshold, mPreStats, mPostStats, iter3->first < 0);
                for(int i=0; i<csPairs.size(); i++) {
                    outputPair(csPairs[i]);
                    delete csPairs[i];
                }
                delete iter3->second;
                iter3 = iter2->second.erase(iter3);
            }
            if(iter2->second.size() == 0) {
                iter2 = iter1->second.erase(iter2);
            } else {
                if(processedPos > iter2->first)
                    processedPos = iter2->first;
                iter2++;
            }
        }
        if(iter1->second.size() == 0) {
            iter1 = mProperClusters.erase(iter1);
            curProcessedPos = processedPos;
        } else {
            if(curProcessedTid > iter1->first) {
                curProcessedTid = iter1->first;
                curProcessedPos = processedPos;
            }
            iter1++;
        }
    }
    if(curProcessedTid != INT_MAX) {
        mProcessedTid = curProcessedTid;
        mProcessedPos = curProcessedPos;
    }
}

void Gencore::finishConsensus(map<int, map<int, unordered_map<long, Cluster*>>>& clusters) {
    map<int, map<int, unordered_map<long, Cluster*>>>::iterator iter1;
    map<int, unordered_map<long, Cluster*>>::iterator iter2;
    unordered_map<long, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end();) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ) {
                if(iter1->first < 0 || iter2->first < 0 ) {
                    unordered_map<string, Pair*>::iterator iterOfPairs;
                    for(iterOfPairs = iter3->second->mPairs.begin(); iterOfPairs!=iter3->second->mPairs.end(); iterOfPairs++) {
                        outputPair(iterOfPairs->second);
                        delete iterOfPairs->second;
                    }
                } else {
                    vector<Pair*> csPairs = iter3->second->clusterByUMI(mOptions->unproperReadsUmiDiffThreshold, mPreStats, mPostStats, iter3->first < 0);
                    for(int i=0; i<csPairs.size(); i++) {
                        outputPair(csPairs[i]);
                        delete csPairs[i];
                    }
                }
                delete iter3->second;
                iter3 = iter2->second.erase(iter3);
            }
            if(iter2->second.size() == 0) {
                iter2 = iter1->second.erase(iter2);
            } else {
                iter2++;
            }
        }
        if(iter1->second.size() == 0) {
            iter1 = clusters.erase(iter1);
        } else {
            iter1++;
        }
    }
}

void Gencore::addToUnProperCluster(bam1_t* b) {
    int tid = b->core.tid;
    int left = b->core.pos;
    long right = b->core.mpos;
    if(b->core.mtid < b->core.tid) {
        tid = b->core.mtid;
        left = b->core.mpos;
        right = b->core.pos;
    }
    createCluster(mUnProperClusters, tid, left, right);
    mUnProperClusters[tid][left][right]->addRead(b);
}

void Gencore::createCluster(map<int, map<int, unordered_map<long, Cluster*>>>& clusters, int tid, int left, long right) {
    auto iter1 = clusters.find(tid);

    if(iter1 == clusters.end()) {
        clusters[tid] = map<int, unordered_map<long, Cluster*>>();
        clusters[tid][left] = unordered_map<long, Cluster*>();
        clusters[tid][left][right] = new Cluster(mOptions);
    } else {
        auto iter2 = iter1->second.find(left);
        if(iter2 == iter1->second.end()) {
            clusters[tid][left] = unordered_map<long, Cluster*>();
            clusters[tid][left][right] = new Cluster(mOptions);
        } else {
            auto iter3 = iter2->second.find(right);
            if(iter3 == iter2->second.end())
                clusters[tid][left][right] = new Cluster(mOptions);
        }
    }
}

void Gencore::addToCluster(bam1_t* b) {
    if(b->core.tid < 0) {
        addToUnProperCluster(b);
    } else {
        addToProperCluster(b);
    }
}
