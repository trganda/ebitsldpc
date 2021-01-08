#ifndef LAB_BINARY_5GLDPC_CODEC_H
#define LAB_BINARY_5GLDPC_CODEC_H

#include "binaryldpccodec.h"
#include "utility.h"
#include <cstdio>
#include <cstring>

namespace lab {
class Binary5GLDPCCodec : public BinaryLDPCCodec {
 public:
  Binary5GLDPCCodec() = default;
  Binary5GLDPCCodec(const Binary5GLDPCCodec &codec);
  explicit Binary5GLDPCCodec(const toml::value &arguments);
  ~Binary5GLDPCCodec() override;
  void Encoder(int *uu, int *cc) const override;
  int Decoder(const double *M2V, int *uu_hat, int iter_cnt) override;

 public:
  int code_len_puncture() const;

 private:
  void SystemMatrixH();

 private:
  int code_len_no_puncture_;// Code length before puncture
  int code_len_puncture_;   // Code length after puncture
  int lifting_factor_;
  int *cc_no_puncture_;
  double *cc_soft_no_puncture_;
};
}// namespace lab
#endif//KMLDPC_BINARY5GLDPCCODEC_H
