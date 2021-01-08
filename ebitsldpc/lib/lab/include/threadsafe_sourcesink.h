#ifndef LAB_THREADSAFE_SOURCESINK_H
#define LAB_THREADSAFE_SOURCESINK_H

#include "sourcesink.h"
#include <condition_variable>
#include <mutex>

namespace lab {
class threadsafe_sourcesink {
 public:
  threadsafe_sourcesink() = default;
  threadsafe_sourcesink(const threadsafe_sourcesink &other) = delete;
  threadsafe_sourcesink &operator=(const threadsafe_sourcesink &other) = delete;
  void GetBitStr(int *uu, int len);
  void ClrCnt();
  void CntErr(const int *uu, const int *uu_hat,
              int len, int accumlator);
  void PrintResult(double snr);

 public:
  std::shared_ptr<unsigned int> try_err_blk();
  std::shared_ptr<unsigned int> try_tot_blk();
  std::shared_ptr<double> try_ber();
  std::shared_ptr<double> try_fer();

 private:
  std::mutex mutex_;
  std::condition_variable data_cond_;
  SourceSink ssink_;
};
}// namespace lab
#endif//LAB_THREADSAFE_SOURCESINK_H
