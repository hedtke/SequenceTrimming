//
//  ConcurrentQueue.h
//  
//
//  Created by Matthias Mueller-Hannemann on 23.07.14.
//
//

#ifndef _ConcurrentQueue_h
#define _ConcurrentQueue_h

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

template <typename T>
class ConcurrentQueue
{
public:
    
    
    bool empty()  {
        std::unique_lock<std::mutex> mlock(mutex_);
        return queue_.empty();
    }
    
    int size()  {
        std::unique_lock<std::mutex> mlock(mutex_);
        return queue_.size();
    }
    
    T pop()
    {
        std::unique_lock<std::mutex> mlock(mutex_);
        while (queue_.empty())
        {
            cond_.wait(mlock);
        }
        auto val = queue_.front();
        queue_.pop();
        return val;
    }
    
    bool tryPop (T& item){
        std::unique_lock<std::mutex> mlock(mutex_);
        if (queue_.empty()) return false;
        
        item = queue_.front();
        queue_.pop();
        
        return true;
    }
    
    void push(const T& item)
    {
        std::unique_lock<std::mutex> mlock(mutex_);
        queue_.push(item);
        mlock.unlock();
        cond_.notify_one();
    }
    
    
    ConcurrentQueue()=default;
    ConcurrentQueue(const ConcurrentQueue&) = delete;   // disable copying
    ConcurrentQueue& operator=(const ConcurrentQueue&) = delete; // disable assignment
    
private:
    std::queue<T> queue_;
    std::mutex mutex_;
    std::condition_variable cond_;
};


#endif
