#include "viterbi.h"

namespace lab {

ViterbiCodec::ViterbiCodec(int constraint, const std::vector<int> &polynomials)
    : constraint_(constraint), polynomials_(polynomials) {
  for (int polynomial : polynomials_) {
    assert(polynomial > 0);
    assert(polynomial < (1 << constraint_));
  }
  initializeOutputs();
}

std::string
ViterbiCodec::Decode(const std::string &bits) const {
  // Compute path metrics and generate trellis.
  Trellis trellis;
  std::vector<int> path_metrics(1 << (constraint_ - 1),
                                std::numeric_limits<int>::max());
  path_metrics.front() = 0;
  for (size_t i = 0; i < bits.size(); i += num_parity_bits()) {
    std::string current_bits(bits, i, num_parity_bits());
    // If some bits are missing, fill with trailing zeros.
    // This is not ideal but it is the best we can do.
    if (current_bits.size() < num_parity_bits()) {
      current_bits.append(
          std::string(num_parity_bits() - current_bits.size(), '0'));
    }
    UpdatePathMetrics(current_bits, path_metrics, trellis);
  }

  // Traceback
  std::string decoded;
  int state = std::distance(path_metrics.begin(), std::min_element(path_metrics.begin(), path_metrics.end()));
  for (int i = trellis.size() - 1; i >= 0; i--) {
    decoded += state >> (constraint_ - 2) ? "1" : "0";
    state = trellis[i][state];
  }
  std::reverse(decoded.begin(), decoded.end());

  // Remove (constraint_ - 1) flushing bits.
  return decoded.substr(0, decoded.size() - constraint_ + 1);
}

void
ViterbiCodec::initializeOutputs() {
  outputs_.resize(1 << constraint_);
  for (int i = 0; i < outputs_.size(); i++) {
    for (int j = 0; j < num_parity_bits(); j++) {
      // Reverse polynomial bits to simplify the calculation, since we use the lsb-current
      int polynomial = ReverseBits(constraint_, polynomials_[j]);
      int input = i;
      int output = 0;
      for (int k = 0; k < constraint_; k++) {
        output ^= (polynomial & 1) & (input & 1);
        polynomial >>= 1;
        input >>= 1;
      }
      outputs_[i] += output ? "1" : "0";
    }
  }
}

void
ViterbiCodec::UpdatePathMetrics(const std::string &bits,
                                std::vector<int> &path_metrics,
                                ViterbiCodec::Trellis &trellis) const {
  std::vector<int> new_path_metrics(path_metrics.size());
  std::vector<int> new_trellis_column(1 << (constraint_ - 1));
  for (size_t i = 0; i < path_metrics.size(); i++) {
    auto p = PathMetric(bits, path_metrics, i);
    new_path_metrics[i] = p.first;
    new_trellis_column[i] = p.second;
  }

  path_metrics = new_path_metrics;
  trellis.push_back(new_trellis_column);
}

std::pair<int, int>
ViterbiCodec::PathMetric(const std::string &bits,
                         const std::vector<int> &prev_path_metrics,
                         int state) const {
  // Find the source state that can transform to state
  int s = (state << 1) & ((1 << (constraint_ - 1)) - 1);
  int source_state1 = s | 0;
  int source_state2 = s | 1;

  int path_metric1 = prev_path_metrics[source_state1];
  if (path_metric1 < std::numeric_limits<int>::max()) {
    path_metric1 += BranchMetric(bits, source_state1, state);
  }
  int path_metric2 = prev_path_metrics[source_state2];
  if (path_metric2 < std::numeric_limits<int>::max()) {
    path_metric2 += BranchMetric(bits, source_state2, state);
  }

  return path_metric1 <= path_metric2 ? std::make_pair(path_metric1, source_state1) :
      std::make_pair(path_metric2, source_state2);
}

int
ViterbiCodec::BranchMetric(
    const std::string &bits,
    int source_state, int target_state) const {
  assert(bits.size() == num_parity_bits());
  assert((target_state & ((1 << (constraint_ - 2)) - 1)) == source_state >> 1);
  const std::string output =
      Output(source_state, target_state >> (constraint_ - 2));

  return HammingDistance(bits, output);
}

std::string
ViterbiCodec::Output(int current_state, int input) const {
  return outputs_.at(current_state | (input << (constraint_ - 1)));
}

int
ReverseBits(int num_bits, int input) {
  int ret = 0;
  while (num_bits-- > 0) {
    ret = (ret << 1) + (input & 1);
    input >>= 1;
  }
  return ret;
}

int
HammingDistance(const std::string &x, const std::string &y) {
  assert(x.size() == y.size());
  int distance = 0;
  for (int i = 0; i < x.size(); i++) {
    distance += x[i] != y[i];
  }
  return distance;
}

}