#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <queue>
#include <string>
#include <utility>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <stdexcept>

using namespace std;

template<typename TTask, typename TResult>
class TaskProcessor {
public:
    virtual ~TaskProcessor() {}
    virtual TResult process(const TTask& task) = 0;
};

template<typename TTask, typename TResult>
class ThreadPool {
public:
    ThreadPool(int threadCount, TaskProcessor<TTask, TResult>* processor);
    ~ThreadPool();

    void submit(const TTask& task);
    TResult take();
    void waitUntilDone();

private:
    TaskProcessor<TTask, TResult>* mProcessor;
    vector<thread> mWorkers;
    queue<TTask> mTaskQueue;
    queue<TResult> mReadyResults;
    mutex mMutex;
    condition_variable mTaskCv;
    condition_variable mResultCv;
    condition_variable mIdleCv;
    bool mStopping;
    bool mFailed;
    string mError;
    size_t mPendingTasks;

    void workerLoop();
};

template<typename TTask, typename TResult>
ThreadPool<TTask, TResult>::ThreadPool(int threadCount, TaskProcessor<TTask, TResult>* processor) {
    mProcessor = processor;
    mStopping = false;
    mFailed = false;
    mPendingTasks = 0;

    for(int i=0; i<threadCount; i++) {
        mWorkers.push_back(thread(&ThreadPool<TTask, TResult>::workerLoop, this));
    }
}

template<typename TTask, typename TResult>
ThreadPool<TTask, TResult>::~ThreadPool() {
    {
        unique_lock<mutex> lock(mMutex);
        mStopping = true;
    }
    mTaskCv.notify_all();
    mResultCv.notify_all();
    mIdleCv.notify_all();
    for(size_t i=0; i<mWorkers.size(); i++) {
        if(mWorkers[i].joinable())
            mWorkers[i].join();
    }
}

template<typename TTask, typename TResult>
void ThreadPool<TTask, TResult>::submit(const TTask& task) {
    unique_lock<mutex> lock(mMutex);
    if(mStopping)
        throw runtime_error("ThreadPool has been stopped");
    if(mFailed)
        throw runtime_error(mError);

    mTaskQueue.push(task);
    mPendingTasks++;
    mTaskCv.notify_one();
}

template<typename TTask, typename TResult>
TResult ThreadPool<TTask, TResult>::take() {
    unique_lock<mutex> lock(mMutex);
    while(true) {
        if(mFailed)
            throw runtime_error(mError);

        if(!mReadyResults.empty()) {
            TResult result = mReadyResults.front();
            mReadyResults.pop();
            return result;
        }

        mResultCv.wait(lock);
    }
}

template<typename TTask, typename TResult>
void ThreadPool<TTask, TResult>::waitUntilDone() {
    unique_lock<mutex> lock(mMutex);
    while(mPendingTasks > 0 && !mFailed) {
        mIdleCv.wait(lock);
    }
    if(mFailed)
        throw runtime_error(mError);
}

template<typename TTask, typename TResult>
void ThreadPool<TTask, TResult>::workerLoop() {
    while(true) {
        TTask taskItem;
        {
            unique_lock<mutex> lock(mMutex);
            while(!mStopping && !mFailed && mTaskQueue.empty()) {
                mTaskCv.wait(lock);
            }
            if(mStopping || mFailed) {
                return;
            }
            taskItem = mTaskQueue.front();
            mTaskQueue.pop();
        }

        TResult result;
        try {
            result = mProcessor->process(taskItem);
        } catch(const exception& e) {
            unique_lock<mutex> lock(mMutex);
            mFailed = true;
            mError = e.what();
            mResultCv.notify_all();
            mIdleCv.notify_all();
            mTaskCv.notify_all();
            return;
        } catch(...) {
            unique_lock<mutex> lock(mMutex);
            mFailed = true;
            mError = "ThreadPool worker failed with unknown exception";
            mResultCv.notify_all();
            mIdleCv.notify_all();
            mTaskCv.notify_all();
            return;
        }

        {
            unique_lock<mutex> lock(mMutex);
            mReadyResults.push(result);
            mPendingTasks--;
            if(mPendingTasks == 0)
                mIdleCv.notify_all();
        }
        mResultCv.notify_all();
    }
}

#endif
