#ifndef LAB_MODEM_H
#define LAB_MODEM_H

#include "log.h"
#include "toml.hpp"
#include "utility.h"
#include <complex>
#include <iostream>
#include <vector>

namespace lab {
class Modem {
 public:
  Modem() = default;
  explicit Modem(const toml::value &arguments);
  virtual ~Modem() = default;

 public:
  void Mapping(const int *bin_cc, std::vector<std::complex<double>> &xx);
  void DeMapping(
      double *bitLin, double *symRin,
      double *bitLout, int yy_len) const;
  std::vector<std::complex<double>> constellations() const;

 private:
  void init(const std::string &modem_file);

 protected:
  // Number of constellation points
  int symbol_num_;
  int input_len_;
  int output_len_;
  std::vector<std::complex<double>> symbol_out_;

 private:
  std::vector<std::vector<int>> symbol_in_;
};
}// namespace lab
#endif
