#ifndef LAB_XOR_SEG_CODEC_H
#define LAB_XOR_SEG_CODEC_H

#include "binary5gldpccodec.h"
#include "binaryldpccodec.h"
#include "log.h"
#include "mat.h"
#include "modemlinearsystem.h"
#include "randnum.h"
#include "utility.h"
#include <algorithm>
#include <complex>
#include <iomanip>
#include <vector>

namespace lab {
class XORSegCodec {
 public:
  XORSegCodec() = default;
  XORSegCodec(const XORSegCodec &codec);
  explicit XORSegCodec(toml::value arguments);
  virtual ~XORSegCodec();
  void Encoder(int *uu, int *cc);
  void Decoder(
      ModemLinearSystem &modem_linear_system,
      const std::vector<std::complex<double>> &hHats, int *uu_hat);
  std::vector<double> GetHistogramData(
      ModemLinearSystem &mlsystem,
      const std::vector<std::complex<double>> &hhats, int *uu_hat);
  int uu_len() const;
  int cc_len() const;

 private:
  void DeMapping(
      ModemLinearSystem &modem_linear_system,
      std::vector<std::pair<int, std::complex<double>>> &thetaList) const;
  int GetParityCheck() const;
  std::vector<double> GetMetrics(
      ModemLinearSystem &modem_linear_system,
      const std::vector<std::complex<double>> &h_hats, int *uu_hat);
  double Metric(int *uu_hat);

 private:
  const toml::value arguments_;
  std::unique_ptr<BinaryLDPCCodec> ldpc_codec_;
  int iter_cnt_;// iteration times while using LDPC on 5G
  bool using_ldpc_5g_;
  bool using_syndrom_metric_;
  int uu_len_;       // length of uu for LDPC
  int cc_len_;       // length of cc for LDPC
  int *rr_;          // hard decision
  double *bit_l_in_; // bit input probability
  double *bit_l_out_;// bit out probability
};
}// namespace lab
#endif