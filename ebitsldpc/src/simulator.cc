#include "simulator.h"

Simulator::Simulator(toml::value &arguments)
    : arguments_(std::move(arguments)), codec_(XORSegCodec(arguments_)),
      codec_data_(codec_.uu_len() + codec_.ebits_len(), codec_.cc_len()),
      modem_linear_system_(lab::ModemLinearSystem(arguments_, codec_.cc_len())) {
  const auto range = toml::find(arguments_, "range");
  min_snr_ = toml::find<double>(range, "minimum_snr");
  max_snr_ = toml::find<double>(range, "maximum_snr");
  step_snr_ = toml::find<double>(range, "step_snr");
  max_err_blk_ = toml::find<int>(range, "maximum_error_number");
  max_num_blk_ = toml::find<int>(range, "maximum_block_number");
  thread_num_blk_ = toml::find<int>(range, "thread_block_number");
  int angle_size = toml::find<int>(range, "angle_size");
  angles_ = std::vector<std::complex<double>>(angle_size);
  for (size_t i = 0; i < angles_.size(); i++) {
    auto angle = ((lab::kPi / 2) / angles_.size()) * i;
    angles_[i] = std::complex<double>(0, exp(angle));
  }
  std::stringstream stream;
  stream << '[' << std::fixed << std::setprecision(3) << min_snr_ << ',' << step_snr_ << ',' << max_snr_ << ']';
  lab::logger::INFO(stream.str(), true);
  stream.str("");
  stream << '[' << "MAX_ERROR_BLK = " << max_err_blk_ << ',' << "MAX_BLK = " << max_num_blk_ << ']';
  lab::logger::INFO(stream.str(), true);
}

void
Simulator::Run() {
  // Threads number
  const auto max_threads = (unsigned long) ((max_snr_ - min_snr_) / step_snr_ + 1);
  // Save simulation results
  std::vector<std::pair<double, double>> ber_result(max_threads);
  std::vector<std::pair<double, double>> fer_result(max_threads);
  std::vector<std::future<std::pair<double, double>>> ber_and_fer(max_threads);
  // Record metric for histogram
  const auto histogram = toml::find(arguments_, "histogram");
  const bool histogram_enable = toml::find<bool>(histogram, "enable");
  lab::ThreadsPool threads_pool(max_threads);
  for (unsigned long i = 0; i < max_threads; i++) {
    ber_and_fer[i] = threads_pool.submit(
        std::bind(
            &Simulator::run, this, std::ref(this->codec_),
            this->modem_linear_system_, std::ref(this->codec_data_),
            (min_snr_ + step_snr_ * i), histogram_enable));
  }
  for (size_t i = 0; i < max_threads; i++) {
    auto result = ber_and_fer[i].get();
    ber_result[i] = std::pair<double, double>((min_snr_ + step_snr_ * i), result.first);
    fer_result[i] = std::pair<double, double>((min_snr_ + step_snr_ * i), result.second);
  }
  std::stringstream stream;
  stream << "BER Result";
  lab::logger::INFO(stream.str(), true);
  stream.str("");
  for (auto item : ber_result) {
    stream << std::fixed << std::setprecision(3) << std::setfill('0') << std::setw(3) << std::right << item.first
           << ' ' << std::setprecision(14) << item.second;
    lab::logger::INFO(stream.str(), true);
    stream.str("");
  }
  stream << "FER Result";
  lab::logger::INFO(stream.str(), true);
  stream.str("");
  for (auto item : fer_result) {
    stream << std::fixed << std::setprecision(3) << std::setfill('0') << std::setw(3) << std::right << item.first
           << ' ' << std::setprecision(14) << item.second;
    lab::logger::INFO(stream.str(), true);
    stream.str("");
  }
}

std::pair<double, double>
Simulator::run(XORSegCodec &codec,
               lab::ModemLinearSystem mls, CodecData &cdata,
               double snr, bool histogram_enable) {
  //var_ = pow(10.0, -0.1 * (snr)) / (codec_.m_coderate * modem_linear_system_.modem_.input_len_);
  double var = pow(10.0, -0.1 * (snr));
  double sigma = sqrt(var);
  mls.set_sigma(sigma);
  mls.set_var(var);
  lab::threadsafe_sourcesink ssink = lab::threadsafe_sourcesink();
  lab::threadsafe_sourcesink ebits_ssink = lab::threadsafe_sourcesink();
  ssink.ClrCnt();
  ebits_ssink.ClrCnt();
  {
    lab::ThreadsPool threads_pool;
    std::vector<std::future<void>> rets;
    auto max_blocks = max_num_blk_;
    auto blocks = 0;
    while (max_blocks > 0) {
      blocks = thread_num_blk_ <= max_blocks ? thread_num_blk_ : max_blocks;
      max_blocks -= blocks;
      rets.push_back(
          threads_pool.submit(
              [this, codec, mls, &ssink, &ebits_ssink, cdata, snr, blocks] {
                run_blocks(
                    codec, mls, std::ref(ssink), std::ref(ebits_ssink), cdata,
                    snr, blocks);
              }));
    }
    // Waiting for the working task to finished
    for (auto &ret : rets) { ret.get(); }
  }
  ssink.PrintResult(snr);
  // BER and FER
  auto ret = std::pair<double, double>(*ssink.try_ber(), *ssink.try_fer());
  return ret;
}

void
Simulator::run_blocks(XORSegCodec codec,
                      lab::ModemLinearSystem mls, lab::threadsafe_sourcesink &ssink,
                      lab::threadsafe_sourcesink &ebits_ssink, CodecData cdata,
                      double snr, const unsigned int max_block) const {
  for (unsigned int i = 0; i < max_block; i++) {
    if (*ssink.try_tot_blk() >= max_num_blk_ || *ssink.try_err_blk() >= max_err_blk_) { return; }
    ssink.GetBitStr(cdata.uu_, cdata.uu_len_);
    ebits_ssink.GetBitStr(cdata.uu_ + codec.uu_len(), codec_.ebits_len());
    codec.Encoder(cdata.uu_, cdata.cc_);

    // Generate H
    std::vector<std::complex<double>> generated_h(codec_.ebits_len());
    for (size_t j = 0; j < generated_h.size(); j++) {
      generated_h[j] = angles_[codec_data_.uu_[codec_data_.uu_len_ + j]];
    }

    // Modulation and pass through the channel
    mls.PartitionModemLSystem(cdata.cc_, generated_h);

    ssink.CntErr(cdata.uu_, cdata.uu_hat_, codec.uu_len(), 1);
    if (int(*ssink.try_tot_blk()) > 0 && int(*ssink.try_tot_blk()) % 100 == 0) { ssink.PrintResult(snr); }
  }
}