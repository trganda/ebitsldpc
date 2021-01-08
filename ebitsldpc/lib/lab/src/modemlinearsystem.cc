#include "modemlinearsystem.h"

namespace lab {
ModemLinearSystem::ModemLinearSystem(const toml::value &arguments, int cc_len)
    : Modem(arguments), cc_len_(cc_len), sigma_(0.0),
      var_(0.0), sym_prob_(nullptr) {
  if (cc_len_ % input_len_ != 0) {
    lab::logger::ERROR(std::string(
        "(cc_len_ = " + std::to_string(cc_len_) + " %% (input_len_ = " + std::to_string((input_len_)) + " ) != 0 !"),
                       true);
    exit(-1);
  }
  xx_ = std::vector<std::complex<double>>(cc_len_ / input_len_);
  yy_ = std::vector<std::complex<double>>(cc_len_ / input_len_);
  sym_prob_ = new double[(cc_len_ / input_len_) * symbol_num_];
}

ModemLinearSystem::ModemLinearSystem(const ModemLinearSystem &mls)
    : Modem(mls), cc_len_(mls.cc_len_),
      sigma_(mls.sigma_), var_(mls.var_),
      xx_(mls.xx_), yy_(mls.yy_) {
  sym_prob_ = new double[(cc_len_ / input_len_) * symbol_num_];
}

ModemLinearSystem::~ModemLinearSystem() {
  delete[] sym_prob_;
}

void
ModemLinearSystem::PartitionModemLSystem(
    const int *cc,
    std::vector<std::complex<double>> &select_h) {
  Mapping(cc, xx_);
  PartitionHAWGNSystem(select_h);
}

void
ModemLinearSystem::PartitionHAWGNSystem(std::vector<std::complex<double>> &h) {
  std::vector<std::complex<double>> noise(xx_.size());
  CLCRandNum::Get().Normal(noise);
  for (size_t i = 0; i < h.size(); i++) {
    auto num_of_part = xx_.size() / h.size();
    for (size_t j = i * num_of_part; j < num_of_part; j++) {
      std::complex<double> temp = xx_[j] * h[i];
      yy_[j] = temp + noise[j] * std::complex<double>(sigma_ / kSqrt2, 0);
    }
  }
}

void
ModemLinearSystem::SoftAWGNDemodulation(
    const std::complex<double> &yy, double *sym_prob,
    std::complex<double> &theta_h) const {
  std::vector<double> symbol_prob(symbol_num_);
  double sqr_norm = 0.0;
  for (size_t i = 0; i < symbol_prob.size(); i++) {
    auto symbol = symbol_out_[i];
    symbol *= theta_h;
    symbol -= yy;
    sqr_norm = (symbol.real() * symbol.real() + symbol.imag() * symbol.imag()) / var_;
    symbol_prob[i] = -sqr_norm;
  }
  auto max_prob = std::max_element(symbol_prob.begin(), symbol_prob.end());
  for (int i = 0; i < symbol_num_; i++) {
    symbol_prob[i] = exp(symbol_prob[i] - *max_prob);
  }
  // normalization
  double sum = 0.0;
  for (int i = 0; i < symbol_num_; i++) {
    sum += symbol_prob[i];
  }
  for (int i = 0; i < symbol_num_; i++) {
    symbol_prob[i] /= sum;
  }
  for (int i = 0; i < symbol_num_; i++) {
    sym_prob[i] = symbol_prob[i];
  }
  utility::ProbClip(sym_prob, symbol_num_);
}

void
ModemLinearSystem::SoftDemodulation(std::vector<std::pair<int, std::complex<double>>> &thetaList) const {
  auto symbolPerPart = yy_.size() / thetaList.size();
  for (size_t i = 0; i < thetaList.size(); i++) {
    for (size_t j = 0; j < symbolPerPart; j++) {
      auto temp = (j + i * symbolPerPart) * symbol_num_;
      SoftAWGNDemodulation(yy_[j + i * symbolPerPart], sym_prob_ + temp, thetaList[i].second);
    }
  }
}

void
ModemLinearSystem::DeMapping(
    std::vector<std::pair<int, std::complex<double>>> &thetaList, double *bitLin,
    double *bitLout) {
  SoftDemodulation(thetaList);
  Modem::DeMapping(bitLin, sym_prob_, bitLout, yy_.size());
}

std::vector<std::complex<double>>
ModemLinearSystem::GetRecvSymbol() const {
  return yy_;
}

void
ModemLinearSystem::set_sigma(double sigma) {
  this->sigma_ = sigma;
}

void
ModemLinearSystem::set_var(double var) {
  this->var_ = var;
}
}// namespace lab