#ifndef LAB_MODEM_LINEAR_SYSTEM_HPP
#define LAB_MODEM_LINEAR_SYSTEM_HPP

#include "modem.h"
#include "randnum.h"
#include "toml.hpp"
#include "utility.h"
#include <complex>
#include <cstring>

namespace lab {
class ModemLinearSystem : public Modem {
 public:
  ModemLinearSystem(const toml::value &arguments, int cc_len);
  ModemLinearSystem(const ModemLinearSystem &mls);
  ~ModemLinearSystem() override;
  void PartitionModemLSystem(
      const int *cc,
      std::vector<std::complex<double>> &select_h);
  void DeMapping(
      std::vector<std::pair<int, std::complex<double>>> &thetaList,
      double *bitLin, double *bitLout);
  std::vector<std::complex<double>> GetRecvSymbol() const;

 public:
  void set_sigma(double sigma);
  void set_var(double var);

 private:
  void SoftDemodulation(std::vector<std::pair<int, std::complex<double>>> &thetaList) const;
  void PartitionHAWGNSystem(std::vector<std::complex<double>> &h);
  void SoftAWGNDemodulation(
      const std::complex<double> &yy, double *sym_prob,
      std::complex<double> &theta_h) const;

 private:
  int cc_len_;
  double sigma_;
  double var_;
  double *sym_prob_;
  std::vector<std::complex<double>> xx_;
  std::vector<std::complex<double>> yy_;
};
}// namespace lab
#endif