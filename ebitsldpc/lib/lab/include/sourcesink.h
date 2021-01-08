#ifndef LAB_SOURCE_SINK_H
#define LAB_SOURCE_SINK_H

#include "log.h"
#include "randnum.h"
#include <complex>
#include <iomanip>
#include <vector>

namespace lab {
class SourceSink {
 public:
  void GetBitStr(int *uu, int len);
  void GetSymStr(int *uu, int qary, int len);
  void ClrCnt();
  void CntErr(
      const int *uu, const int *uu_hat,
      int len, int accumulator);
  void PrintResult(double snr) const;
  // Getter
  unsigned int err_blk() const;
  unsigned int tot_blk() const;
  double ber() const;
  double fer() const;

 private:
  unsigned int tot_blk_;
  unsigned int tot_bit_;
  unsigned int err_blk_;
  unsigned int err_bit_;
  double ber_;
  double fer_;
};
}// namespace lab
#endif