#include "sourcesink.h"

namespace lab {
void
SourceSink::GetBitStr(int *uu, int len) {
  for (int t = 0; t < len; t++) {
    uu[t] = (CLCRandNum::Get().Uniform() < 0.5 ? 0 : 1);
  }
}

void
SourceSink::GetSymStr(int *uu, int qary, int len) {
  for (int t = 0; t < len; t++) {
    uu[t] = qary;
    while (uu[t] == qary)
      uu[t] = (int) (qary * CLCRandNum::Get().Uniform());
  }
}

void
SourceSink::ClrCnt() {
  tot_blk_ = 0;
  tot_bit_ = 0;
  err_blk_ = 0;
  err_bit_ = 0;
}

void
SourceSink::CntErr(
    const int *uu, const int *uu_hat,
    int len, int accumulator) {
  int temp_err = 0;
  for (int t = 0; t < len; t++) {
    if (uu_hat[t] != uu[t])
      temp_err++;
  }
  if (accumulator == 1) {
    if (temp_err > 0) {
      err_bit_ += temp_err;
      err_blk_ += 1;
    }
    tot_blk_ += 1.0;
    tot_bit_ += len;
    ber_ = double(err_bit_) / tot_bit_;
    fer_ = double(err_blk_) / tot_blk_;
  }
}

void
SourceSink::PrintResult(double snr) const {
  std::stringstream stream;
  stream << std::fixed << std::setprecision(3) << std::setfill('0')
         << "SNR = "
         << std::setw(3) << std::right << snr << ' '
         << "Total blk = "
         << std::setw(7) << std::right << std::setprecision(0) << tot_blk_ << ' '
         << "Error blk = "
         << std::setw(7) << std::right << err_blk_ << ' '
         << "Error bit = "
         << std::setw(7) << std::right << err_bit_ << ' '
         << std::fixed << std::setprecision(14)
         << "BER = " << ber_ << ' '
         << "FER = " << fer_;
  logger::INFO(stream.str(), true);
}

unsigned int
SourceSink::tot_blk() const {
  return tot_blk_;
}

unsigned int
SourceSink::err_blk() const {
  return err_blk_;
}

double
SourceSink::ber() const {
  return ber_;
}

double
SourceSink::fer() const {
  return fer_;
}
}// namespace lab
