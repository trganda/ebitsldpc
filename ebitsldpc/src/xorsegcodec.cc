#include "xorsegcodec.h"

XORSegCodec::XORSegCodec(const XORSegCodec &other)
    : arguments_(other.arguments_), iter_cnt_(other.iter_cnt_),
      using_ldpc_5g_(other.using_ldpc_5g_), using_syndrom_metric_(other.using_syndrom_metric_),
      uu_len_(other.uu_len_), ebits_len_(other.ebits_len_), cc_len_(other.cc_len_),
      generator_matrix_(other.generator_matrix_) {
  if (using_ldpc_5g_) {
    // For unknown reason, the copy constructor of Binary5GLDPCCodec will
    // occurs some memory problem. It's recommend to  use the one-arguments
    // constructor to build a new object.
    ldpc_codec_ = std::make_unique<lab::Binary5GLDPCCodec>(arguments_);
  } else {
    ldpc_codec_ = std::make_unique<lab::BinaryLDPCCodec>(*other.ldpc_codec_);
  }
  rr_ = new int[cc_len_];
  bit_l_in_ = new double[cc_len_];
  bit_l_out_ = new double[cc_len_];
}

XORSegCodec::XORSegCodec(toml::value arguments)
    : arguments_(std::move(arguments)) {
  const auto xcodec = toml::find(arguments_, "xcodec");
  using_ldpc_5g_ = toml::find<bool>(xcodec, "5gldpc");
  using_syndrom_metric_ = toml::find<bool>(xcodec, "metric_type");
  iter_cnt_ = toml::find<int>(xcodec, "metric_iter");
  ebits_len_ = toml::find<int>(xcodec, "ebits_length");
  if (using_ldpc_5g_) {
    lab::logger::INFO("Using 5G LDPC.", true);
    auto temp = new lab::Binary5GLDPCCodec(arguments_);
    cc_len_ = temp->code_len_puncture();
    ldpc_codec_ = std::unique_ptr<lab::Binary5GLDPCCodec>(temp);
  } else {
    lab::logger::INFO("Using traditional LDPC.", true);
    ldpc_codec_ = std::make_unique<lab::BinaryLDPCCodec>(arguments_);
    cc_len_ = ldpc_codec_->code_len();
  }
  // Random generator matrix
  generator_matrix_ = std::vector<std::vector<int>>(ebits_len_, std::vector<int>(cc_len_, 0));
  for (auto &column : generator_matrix_) {
    for (auto &item : column) {
      item = lab::CLCRandNum::Get().Uniform() < 0.5 ? 0 : 1;
    }
  }
  uu_len_ = ldpc_codec_->code_dim();
  rr_ = new int[cc_len_];
  bit_l_in_ = new double[cc_len_];
  bit_l_out_ = new double[cc_len_];
}

XORSegCodec::~XORSegCodec() {
  delete[] rr_;
  delete[] bit_l_in_;
  delete[] bit_l_out_;
}

void
XORSegCodec::Encoder(int *uu, int *cc) {
  ldpc_codec_->Encoder(uu, cc);
  int temp[cc_len_];
  lab::utility::MatrixProd(uu + uu_len_, temp,
                           generator_matrix_, ebits_len_, cc_len_);
  for (int i = 0; i < cc_len_; i++) {
    cc[i] ^= temp[i];
  }
}

void
XORSegCodec::Decoder(lab::ModemLinearSystem &mls, std::vector<std::complex<double>> &angles, int *uu_hat) {
  DecodeEbits(mls, angles, uu_hat);
  int temp[cc_len_];
  lab::utility::MatrixProd(uu_hat + uu_len_, temp,
                           generator_matrix_, ebits_len_, cc_len_);
  for (int i=0; i<cc_len_; i++) {
    if (temp[i]) {
      bit_l_out_[i] = 1 - bit_l_out_[i];
    }
  }
  ldpc_codec_->Decoder(bit_l_out_, uu_hat, ldpc_codec_->max_iter());
}

int
XORSegCodec::uu_len() const {
  return uu_len_;
}

int
XORSegCodec::ebits_len() const {
  return ebits_len_;
}

int
XORSegCodec::cc_len() const {
  return cc_len_;
}

void
XORSegCodec::DecodeEbits(lab::ModemLinearSystem &mls, std::vector<std::complex<double>> &angles, int *uu_hat) {
  std::vector<std::vector<double>> branch_metrics(ebits_len_, std::vector<double>(angles.size()));
  auto received_symbols = mls.GetRecvSymbol();
  auto constellations = mls.constellations();
  auto partition_length = received_symbols.size() / ebits_len_;
  for (size_t i=0; i<branch_metrics.size(); i++) {
    auto symbols = std::vector<std::complex<double>>(received_symbols.begin() + i * partition_length,
        received_symbols.begin() + (i+1) * partition_length);
    for (size_t j=0; j<branch_metrics[i].size(); j++) {
      branch_metrics[i][j] = AdaptFunc(symbols, constellations, angles[j], mls.var());
    }
  }
}

void
XORSegCodec::DeMapping(
    lab::ModemLinearSystem &mls,
    std::vector<std::pair<int, std::complex<double>>> &thetaList) const {
  // demapping to Get soft information
  for (int i = 0; i < cc_len_; i++) {
    bit_l_in_[i] = 0.5;
  }
  mls.DeMapping(
      thetaList,
      bit_l_in_, bit_l_out_);
}

double
XORSegCodec::AdaptFunc(std::vector<std::complex<double>> &data,
                       std::vector<std::complex<double>> &constellations,
                       std::complex<double> &h, double var) {
  double result = 0.0;
  double co = 1.0 / sqrt(2 * lab::kPi * var);

  auto symbols = constellations;
  for (auto &symbol : symbols) {
    symbol *= h;
  }

  for (auto i : data) {
    double sum = 0.0;
    for (auto j : symbols) {
      double alpha = 1.0 / symbols.size();
      std::complex<double> differ = i - j;
      sum += alpha * co * exp(-((pow(differ.real(), 2) + pow(differ.imag(), 2)) / var));
    }
    result += log(sum);
  }
  result /= double(data.size());
  return result;
}