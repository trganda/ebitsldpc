#include "modem.h"

namespace lab {
Modem::Modem(const toml::value &arguments)
    : symbol_num_(0), input_len_(0), output_len_(0) {
  const auto modem = toml::find(arguments, "modem");
  std::string modem_file = toml::find<std::string>(modem, "modem_file");
  this->init(modem_file);
}

void
Modem::Mapping(const int *bin_cc, std::vector<std::complex<double>> &xx) {
  for (size_t i = 0; i < xx.size(); i++) {
    int symbol_idx = 0;
    for (int j = 0; j < input_len_; j++) {
      symbol_idx = (symbol_idx << 1) + bin_cc[j + i * input_len_];
    }
    xx[i] = symbol_out_[symbol_idx];
  }
}

void
Modem::DeMapping(
    double *bitLin, double *symRin,
    double *bitLout, int yy_len) const {
  utility::ProbClip(bitLin, yy_len * input_len_);
  utility::ProbClip(symRin, symbol_num_ * yy_len);
  std::vector<double> symbol_prob(symbol_num_);
  int t = 0;
  for (int i = 0; i < yy_len; i++) {
    for (double &j : symbol_prob) {
      j = 1.0;
    }
    // Convert bit extrinsic message to symbol message
    t = i * input_len_;
    for (int j = 0; j < input_len_; j++) {
      for (int k = 0; k < symbol_num_; k++) {
        if (symbol_in_[k][j] == 0) {
          symbol_prob[k] *= bitLin[t];
        } else {
          symbol_prob[k] *= 1.0 - bitLin[t];
        }
      }
      t++;
    }
    // Compute symbol full message
    double sum = 0.0;
    t = i * symbol_num_;
    for (double &j : symbol_prob) {
      j *= symRin[t];
      sum += j;
      t++;
    }
    // Normalized
    for (double &j : symbol_prob) {
      j /= sum;
    }
    // Compute bit extrinsic message for output
    t = i * input_len_;
    double prob_0;
    double prob_1;
    for (int j = 0; j < input_len_; j++) {
      prob_0 = 0.0;
      prob_1 = 0.0;
      for (int k = 0; k < symbol_num_; k++) {
        if (symbol_in_[k][j] == 0) {
          prob_0 += symbol_prob[k];
        } else {
          prob_1 += symbol_prob[k];
        }
      }
      prob_0 /= bitLin[t];
      prob_1 /= (1.0 - bitLin[t]);
      bitLout[t] = prob_0 / (prob_0 + prob_1);
      t++;
    }
  }
  utility::ProbClip(bitLout, yy_len * input_len_);
}

std::vector<std::complex<double>>
Modem::constellations() const {
  return symbol_out_;
}

void
Modem::init(const std::string &modem_file) {
  std::ifstream ifs(modem_file, std::ios_base::binary);
  if (!ifs.is_open()) {
    lab::logger::ERROR("Cannot Opne" + modem_file, true);
    exit(-1);
  }
  std::string temp;
  ifs >> temp;
  ifs >> input_len_;
  ifs >> temp;
  ifs >> output_len_;
  ifs >> temp;
  symbol_num_ = 1 << input_len_;
  symbol_in_ = std::vector<std::vector<int>>(
      symbol_num_,
      std::vector<int>(input_len_));
  symbol_out_ = std::vector<std::complex<double>>(symbol_num_);
  int sym_dec = 0;
  double energies = 0;
  for (int i = 0; i < symbol_num_; i++) {
    ifs >> sym_dec;
    int temp_dec = 0;
    for (int j = 0; j < input_len_; j++) {
      ifs >> symbol_in_[i][j];
      temp_dec = (temp_dec << 1) + symbol_in_[i][j];
    }
    if (sym_dec != temp_dec || sym_dec != i) {
      lab::logger::ERROR(std::string(
          std::to_string(sym_dec) + " is not the binary expression of " + std::to_string(temp_dec)),
                         true);
      exit(-1);
    }
    double real, imag;
    ifs >> real >> imag;
    symbol_out_[i] = std::complex<double>(real, imag);
    energies += pow(abs(symbol_out_[i]), 2);
  }
  ifs.close();
  energies /= symbol_num_;
  for (int i = 0; i < symbol_num_; i++) {
    symbol_out_[i] /= sqrt(energies);
  }
}
}// namespace lab