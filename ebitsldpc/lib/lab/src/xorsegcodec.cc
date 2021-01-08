#include "xorsegcodec.h"

namespace lab {
XORSegCodec::XORSegCodec(const XORSegCodec &codec)
    : arguments_(codec.arguments_), iter_cnt_(codec.iter_cnt_), using_ldpc_5g_(codec.using_ldpc_5g_),
      using_syndrom_metric_(codec.using_syndrom_metric_),
      uu_len_(codec.uu_len_), cc_len_(codec.cc_len_) {
  if (using_ldpc_5g_) {
    // For unknown reason, the copy constructor of Binary5GLDPCCodec will
    // occurs some memory problem. It's recommend to  use the one-arguments
    // constructor to build a new object.
    ldpc_codec_ = std::make_unique<Binary5GLDPCCodec>(arguments_);
  } else {
    ldpc_codec_ = std::make_unique<BinaryLDPCCodec>(*codec.ldpc_codec_);
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
  if (using_ldpc_5g_) {
    logger::INFO("Using 5G LDPC.", true);
    auto temp = new Binary5GLDPCCodec(arguments_);
    cc_len_ = temp->code_len_puncture();
    ldpc_codec_ = std::unique_ptr<Binary5GLDPCCodec>(temp);
  } else {
    logger::INFO("Using traditional LDPC.", true);
    ldpc_codec_ = std::make_unique<BinaryLDPCCodec>(arguments_);
    cc_len_ = ldpc_codec_->code_len();
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
}

void
XORSegCodec::Decoder(
    ModemLinearSystem &modem_linear_system,
    const std::vector<std::complex<double>> &hHats, int *uu_hat) {
  std::vector<double> metric_results;
  std::vector<std::pair<int, std::complex<double>>> temp;
  metric_results = GetMetrics(modem_linear_system, hHats, uu_hat);
  auto minIndex = std::distance(
      metric_results.begin(),
      min_element(metric_results.begin(), metric_results.end()));
  logger::INFO(std::string("hatIndex = " + std::to_string(minIndex)), false);
  temp = {std::pair<int, std::complex<double>>(0, hHats[minIndex])};
  DeMapping(modem_linear_system, temp);
  ldpc_codec_->Decoder(bit_l_out_, uu_hat, ldpc_codec_->max_iter());
}

std::vector<double>
XORSegCodec::GetHistogramData(
    ModemLinearSystem &mlsystem,
    const std::vector<std::complex<double>> &hhats, int *uu_hat) {
  return GetMetrics(mlsystem, hhats, uu_hat);
}

int
XORSegCodec::uu_len() const {
  return uu_len_;
}

int
XORSegCodec::cc_len() const {
  return cc_len_;
}

void
XORSegCodec::DeMapping(
    ModemLinearSystem &modem_linear_system,
    std::vector<std::pair<int, std::complex<double>>> &thetaList) const {
  // demapping to Get soft information
  for (int i = 0; i < cc_len_; i++) {
    bit_l_in_[i] = 0.5;
  }
  modem_linear_system.DeMapping(
      thetaList,
      bit_l_in_, bit_l_out_);
}

int
XORSegCodec::GetParityCheck() const {
  if (using_ldpc_5g_) {
    return ldpc_codec_->ParityCheck(ldpc_codec_->cc_hat());
  } else {
    // hard decision without decoding
    for (int i = 0; i < cc_len_; i++) {
      if (bit_l_out_[i] > 0.5) {
        rr_[i] = 1;
      } else {
        rr_[i] = 0;
      }
    }
    return ldpc_codec_->ParityCheck(rr_);
  }
}

std::vector<double>
XORSegCodec::GetMetrics(
    ModemLinearSystem &modem_linear_system,
    const std::vector<std::complex<double>> &h_hats, int *uu_hat) {
  std::vector<double> metric_results(h_hats.size(), 0);
  std::vector<std::pair<int, std::complex<double>>> temp;
  std::stringstream stream;
  for (size_t i = 0; i < metric_results.size(); i++) {
    temp = {std::pair<int, std::complex<double>>(0, h_hats[i])};
    DeMapping(modem_linear_system, temp);
    metric_results[i] = Metric(uu_hat);
    stream << std::fixed << std::setprecision(14)
           << "Hhat = " << h_hats[i]
           << " Metric = "
           << std::setw(5) << std::right
           << metric_results[i];
    logger::INFO(stream.str(), false);
    stream.str("");
    metric_results[i] = abs(metric_results[i]);
  }
  return metric_results;
}

double
XORSegCodec::Metric(int *uu_hat) {
  double metric_result = 0.0;
  if (using_syndrom_metric_ == 1) {
    ldpc_codec_->Decoder(bit_l_out_, uu_hat, iter_cnt_);
    std::vector<double> soft_syndroms;
    soft_syndroms = std::vector<double>(
        ldpc_codec_->syndrom_soft(),
        ldpc_codec_->syndrom_soft() + ldpc_codec_->num_row());
    for (auto j = 0; j < ldpc_codec_->num_row(); j++) {
      metric_result += log(soft_syndroms[j]);
    }
  } else {
    if (using_ldpc_5g_) {
      ldpc_codec_->Decoder(bit_l_out_, uu_hat, iter_cnt_);
    }
    metric_result = GetParityCheck();
  }
  return metric_result;
}
}// namespace lab