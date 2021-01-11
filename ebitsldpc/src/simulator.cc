#include "simulator.h"

Simulator::Simulator(toml::value arguments)
    : arguments_(std::move(arguments)), codec_(lab::XORSegCodec(arguments_)),
      codec_data_(codec_.uu_len(), codec_.cc_len()),
      modem_linear_system_(lab::ModemLinearSystem(arguments_, codec_.cc_len())) {
  const auto range = toml::find(arguments_, "range");
  min_snr_ = toml::find<double>(range, "minimum_snr");
  max_snr_ = toml::find<double>(range, "maximum_snr");
  step_snr_ = toml::find<double>(range, "step_snr");
  max_err_blk_ = toml::find<int>(range, "maximum_error_number");
  max_num_blk_ = toml::find<int>(range, "maximum_block_number");
  thread_num_blk_ = toml::find<int>(range, "thread_block_number");
  const auto decoder = toml::find(arguments_, "decoder");
  known_h_ = toml::find<bool>(decoder, "true_h_arg");
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
Simulator::run(
    lab::XORSegCodec &codec, lab::ModemLinearSystem mls, CodecData &cdata,
    double snr, bool histogram_enable) {
  //var_ = pow(10.0, -0.1 * (snr)) / (codec_.m_coderate * modem_linear_system_.modem_.input_len_);
  double var = pow(10.0, -0.1 * (snr));
  double sigma = sqrt(var);
  mls.set_sigma(sigma);
  mls.set_var(var);
  lab::threadsafe_sourcesink ssink = lab::threadsafe_sourcesink();
  ssink.ClrCnt();
  std::fstream out;
  if (histogram_enable) {
    std::string histfilename = "histogram_" + std::to_string(snr) + ".txt";
    out = std::fstream(histfilename, std::ios::out);
  }
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
              [this, codec, mls, &ssink, cdata, &out, snr, histogram_enable, blocks] {
                run_blocks(
                    codec, mls, std::ref(ssink), cdata,
                    std::ref(out), snr, histogram_enable, blocks);
              }));
    }
    // Waiting for the working task to finished
    for (auto &ret : rets) { ret.get(); }
  }
  if (histogram_enable) { out.close(); }
  ssink.PrintResult(snr);
  // BER and FER
  auto ret = std::pair<double, double>(*ssink.try_ber(), *ssink.try_fer());
  return ret;
}

void
Simulator::run_blocks(
    lab::XORSegCodec codec, lab::ModemLinearSystem mls, lab::threadsafe_sourcesink &ssink,
    CodecData cdata, std::fstream &out, double snr, bool histogram_enable,
    const unsigned int max_block) const {
  for (unsigned int i = 0; i < max_block; i++) {
    if (*ssink.try_tot_blk() >= max_num_blk_ || *ssink.try_err_blk() >= max_err_blk_) { return; }
    ssink.GetBitStr(cdata.uu_, cdata.uu_len_);
    codec.Encoder(cdata.uu_, cdata.cc_);
    // Generate H
    std::complex<double> true_h;
    lab::CLCRandNum::Get().Normal(true_h);
    true_h *= sqrt(0.5);
    std::stringstream stream;
    stream << "Generated H = " << true_h;
    lab::logger::INFO(stream.str(), false);
    std::vector<std::complex<double>> generated_h(1);
    for (auto &item : generated_h) { item = true_h; }
    // Modulation and pass through the channel
    mls.PartitionModemLSystem(cdata.cc_, generated_h);
    // Get constellation
    auto constellations = mls.constellations();
    // Get received symbols
    auto received_symbols = mls.GetRecvSymbol();
    std::vector<std::complex<double>> h_hats;
    if (known_h_) {
      h_hats.push_back(true_h);
    } else {

    }
    if (histogram_enable) {
      auto metrics = codec.GetHistogramData(mls, h_hats, cdata.uu_hat_);
      auto idx_of_min = std::distance(metrics.begin(), min_element(metrics.begin(), metrics.end()));
      std::mutex mutex;
      std::lock_guard<std::mutex> lock(mutex);
      for (size_t j = idx_of_min; j < idx_of_min + metrics.size(); j++) {
        out << metrics[j % metrics.size()] << ' ';
      }
      out << std::endl;
    } else {
      codec.Decoder(mls, h_hats, cdata.uu_hat_);
    }
    ssink.CntErr(cdata.uu_, cdata.uu_hat_, codec.uu_len(), 1);
    if (int(*ssink.try_tot_blk()) > 0 && int(*ssink.try_tot_blk()) % 100 == 0) { ssink.PrintResult(snr); }
  }
}