#ifndef EBITSLDPC_XORCODEC_H_
#define EBITSLDPC_XORCODEC_H_

#include <algorithm>
#include <complex>
#include <iomanip>
#include <vector>

#include "binary5gldpccodec.h"
#include "binaryldpccodec.h"
#include "log.h"
#include "mat.h"
#include "modemlinearsystem.h"
#include "randnum.h"
#include "utility.h"

class XORSegCodec {
 public:
  XORSegCodec() = default;
  XORSegCodec(const XORSegCodec &other);
  explicit XORSegCodec(toml::value arguments);
  virtual ~XORSegCodec();
  void Encoder(int *uu, int *cc);
  void Decoder(lab::ModemLinearSystem &mls,
               std::vector<std::complex<double>> &angles, int *uu_hat);
  int uu_len() const;
  int ebits_len() const;
  int cc_len() const;

 private:
  void DecodeEbits(lab::ModemLinearSystem &mls, std::vector<std::complex<double>> &angles, int *uu_hat);
  void DeMapping(lab::ModemLinearSystem &mls,
                 std::vector<std::pair<int, std::complex<double>>> &thetaList) const;
  int GetParityCheck() const;
  double AdaptFunc(std::vector<std::complex<double>> &data,
                   std::vector<std::complex<double>> &constellations,
                   std::complex<double> &h,
                   double var);

 private:
  const toml::value arguments_;
  std::unique_ptr<lab::BinaryLDPCCodec> ldpc_codec_;
  int iter_cnt_;    // iteration times while using LDPC on 5G
  bool using_ldpc_5g_;
  bool using_syndrom_metric_;
  int uu_len_;      // length of uu for LDPC
  int ebits_len_;   // length of ebits
  int cc_len_;      // length of cc for LDPC
  int *rr_;         // hard decision
  std::vector<std::vector<int>> generator_matrix_;
  double *bit_l_in_; // bit input probability
  double *bit_l_out_;// bit out probability
};

#endif//EBITSLDPC_XORCODEC_H_
