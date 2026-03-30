#include "reference.h"
#include "util.h"
#include <atomic>

Reference* Reference::mInstance = NULL;

Reference* Reference::instance(Options* opt) {
    if(mInstance == NULL)
        mInstance = new Reference(opt);
    
    return mInstance;
}

Reference::Reference(Options* opt) {
    mOptions = opt;
    mRef = NULL;
    if(!mOptions->refFile.empty()) {
        mRef = new FastaReader(mOptions, mOptions->refFile);
        mRef->readAll();
    }
}

Reference::~Reference() {
    if(mRef) {
        delete mRef;
        mRef = NULL;
    }
    mInstance = NULL;
}

const unsigned char* Reference::getData(int bamContig, int pos, int len) {
    if(mRef == NULL)
        return NULL;
    if(mOptions->bamHeader == NULL)
        return NULL;

    // get contig name from bam header
    string contigName(mOptions->bamHeader->target_name[bamContig]);
    map<string, const unsigned char*>::const_iterator contigIter = mRef->mAllContigs.find(contigName);
    map<string, long>::const_iterator sizeIter = mRef->mAllContigSizes.find(contigName);

    if(contigIter == mRef->mAllContigs.end() || sizeIter == mRef->mAllContigSizes.end()) {
        static atomic<bool> reported(false);
        if(!reported.exchange(true))
            cerr << "contig " << contigName << " not found in the reference, please make sure your reference is correct" << endl;
        return NULL;
    }

    if(pos + len >= sizeIter->second){
        static atomic<bool> reported(false);
        if(!reported.exchange(true))
            cerr << "contig " << contigName << " doesn't match the length in the reference, please make sure your reference is correct" << endl;
        return NULL;
    }

    return contigIter->second;
}
