#ifndef LAB_BINARY_LDPC_CODEC_H
#define LAB_BINARY_LDPC_CODEC_H

#include "log.h"
#include "toml.hpp"
#include "utility.h"
#include <cstdio>
#include <cstring>

namespace lab {
class BinaryLDPCCodec {
 public:
  BinaryLDPCCodec() = default;
  BinaryLDPCCodec(const BinaryLDPCCodec &codec);
  explicit BinaryLDPCCodec(const toml::value &arguments);
  virtual ~BinaryLDPCCodec();

 public:
  virtual void Encoder(int *uu, int *cc) const;
  virtual int Decoder(const double *M2V, int *uu_hat, int iter_count);
  int ParityCheck(const int *rr) const;
  void InitMsg() const;

 public:
  int code_dim() const;
  int code_len() const;
  int num_row() const;
  double *syndrom_soft() const;
  int *cc_hat() const;
  int max_iter() const;

 private:
  virtual void SystemMatrixH();

 protected:
  void FreeTannerGraph() const;

 protected:
  // Code parameters
  int code_dim_;       // Dimension of code
  int code_len_;       // Length of code
  int code_chk_;       // Number of Parity check
  double coderate_;    // Code rate
  bool encoder_active_;// 0： not use encoder  1： otherwise
  // Parity-check matrix
  int num_row_; // Row of parity matrix
  int num_col_; // Column of parity matrix
  char **dec_h_;// Matrix used for decoding
  char **enc_h_;// Matrix used for encoding
  double *syndrom_soft_;
  Edge *row_head_;
  Edge *col_head_;
  int *cc_hat_;
  int max_iter_;
  int success_;
};
}// namespace lab
#endif