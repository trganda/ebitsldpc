#ifndef LAB_THREADPOOL_H_
#define LAB_THREADPOOL_H_

#include "join_threads.h"
#include "threadsafe_queue.h"
#include <atomic>
#include <functional>
#include <future>
#include <thread>
#include <type_traits>

namespace lab {
class ThreadsPool {
 public:
  ThreadsPool() : done_(false), joiner_(threads_) {
    unsigned const thread_count = std::thread::hardware_concurrency();
    try {
      for (unsigned i = 0; i < thread_count; ++i) {
        threads_.emplace_back(&ThreadsPool::worker_thread, this);
      }
    } catch (...) {
      done_ = true;
      throw;
    }
  }

  explicit ThreadsPool(unsigned const thread_count) : done_(false), joiner_(threads_) {
    try {
      for (unsigned i = 0; i < thread_count; ++i) {
        threads_.emplace_back(&ThreadsPool::worker_thread, this);
      }
    } catch (...) {
      done_ = true;
      throw;
    }
  }

  ~ThreadsPool() {
    done_ = true;
  }

  template<typename FunctionType, typename... Args>
  std::future<typename std::result_of<FunctionType(Args...)>::type> submit(FunctionType &&f, Args &&...args) {
    using return_type = typename std::result_of<FunctionType(Args...)>::type;

    auto task = std::make_shared<std::packaged_task<return_type()>>(
        std::bind(std::forward<FunctionType>(f), std::forward<Args>(args)...));
    std::future<return_type> res = task->get_future();
    if (done_) {
      throw std::runtime_error("Push task on stopped thread pool");
    }
    work_queue_.push(
        [task]() {
          (*task)();
        });
    return res;
  }

  void run_pending_task() {
    std::function<void()> task;
    if (work_queue_.try_pop(task)) {
      task();
    } else {
      std::this_thread::yield();
    }
  }

 private:
  void worker_thread() {
    while (!done_) {
      run_pending_task();
    }
  }

 private:
  // Do not change the order of the declared variable
  std::atomic_bool done_;
  threadsafe_queue<std::function<void()>> work_queue_;
  std::vector<std::thread> threads_;
  JoinThreads joiner_;
};
}// namespace lab

#endif
