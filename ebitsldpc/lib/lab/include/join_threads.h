#ifndef LAB_JOIN_THREADS_H_
#define LAB_JOIN_THREADS_H_

#include <thread>
#include <vector>

namespace lab {
class JoinThreads {
 public:
  explicit JoinThreads(std::vector<std::thread> &threads)
      : threads_(threads) {}

  ~JoinThreads() {
    for (auto &thread : threads_) {
      if (thread.joinable())
        thread.join();
    }
  }

 private:
  std::vector<std::thread> &threads_;
};
}// namespace lab
#endif
