#include "threadsafe_sourcesink.h"

namespace lab {
void
threadsafe_sourcesink::GetBitStr(int *uu, int len) {
  std::lock_guard<std::mutex> lock(mutex_);
  ssink_.GetBitStr(uu, len);
}

void
threadsafe_sourcesink::ClrCnt() {
  std::lock_guard<std::mutex> lock(mutex_);
  ssink_.ClrCnt();
}

void
threadsafe_sourcesink::CntErr(const int *uu, const int *uu_hat,
                              int len, int accumlator) {
  std::lock_guard<std::mutex> lock(mutex_);
  ssink_.CntErr(uu, uu_hat, len, accumlator);
}

void
threadsafe_sourcesink::PrintResult(double snr) {
  std::lock_guard<std::mutex> lock(mutex_);
  ssink_.PrintResult(snr);
}

std::shared_ptr<unsigned int>
threadsafe_sourcesink::try_err_blk() {
  std::lock_guard<std::mutex> lock(mutex_);
  return std::make_shared<unsigned int>(ssink_.err_blk());
}

std::shared_ptr<unsigned int>
threadsafe_sourcesink::try_tot_blk() {
  std::lock_guard<std::mutex> lock(mutex_);
  return std::make_shared<unsigned int>(ssink_.tot_blk());
}

std::shared_ptr<double>
threadsafe_sourcesink::try_ber() {
  std::lock_guard<std::mutex> lock(mutex_);
  return std::make_shared<double>(ssink_.ber());
}

std::shared_ptr<double>
threadsafe_sourcesink::try_fer() {
  std::lock_guard<std::mutex> lock(mutex_);
  return std::make_shared<double>(ssink_.fer());
}
}// namespace lab