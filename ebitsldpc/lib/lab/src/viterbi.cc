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

std::vector<std::string>
ViterbiCodec::PLVDecode(const std::string &bits, unsigned int L) const {
  PLV_Trellis trellis;
  // k_path_metrics[j][k] means that k-th path metric to state j
  std::vector<std::vector<int>> k_path_metrics(
      1 << (constraint_ - 1), std::vector<int>(L, std::numeric_limits<int>::max()));
  k_path_metrics.front().front() = 0;

  for (size_t i = 0; i < bits.size(); i += num_parity_bits()) {
    std::string current_bits(bits, i, num_parity_bits());
    // If some bits are missing, fill with trailing zeros.
    // This is not ideal but it is the best we can do.
    if (current_bits.size() < num_parity_bits()) {
      current_bits.append(
          std::string(num_parity_bits() - current_bits.size(), '0'));
    }
    UpdatePathMetrics(current_bits, k_path_metrics, trellis);
  }

  // Traceback
  std::vector<std::string> decodes;
  std::vector<int> end_state(L, 0);
  unsigned int kL = 0;
  for (auto &k_path_metric : k_path_metrics) {
    if (k_path_metric[kL] == std::numeric_limits<int>::max()) {
      break;
    }
    kL++;
  }
  // Temp metrics for finding the last state occupied by k-th best path.
  //  std::vector<int> metrics(k_path_metrics.size());
  //  for (unsigned int i = 0; i < L; i++) {
  //    for (size_t j=0; j<k_path_metrics.size(); j++) {
  //      metrics[j] = k_path_metrics[j][i];
  //    }
  //    end_state[i] = std::distance(metrics.begin(), std::min_element(metrics.begin(), metrics.end()));
  //  }

  for (unsigned int i = 0; i < kL; i++) {
    std::string k_result;
    for (int j = trellis.size() - 1; j >= 0; j--) {
      k_result += end_state[i] >> (constraint_ - 2) ? "1" : "0";
      end_state[i] = trellis[j][end_state[i]][i];
    }
    std::reverse(k_result.begin(), k_result.end());
    decodes.push_back(k_result.substr(0, k_result.size() - constraint_ + 1));
  }
  return decodes;
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
ViterbiCodec::UpdatePathMetrics(
    const std::string &bits,
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

void
ViterbiCodec::UpdatePathMetrics(const std::string &bits,
                                std::vector<std::vector<int>> &k_path_metrics,
                                PLV_Trellis &trellis) const {
  std::vector<std::vector<int>> new_k_path_metrics(
      k_path_metrics.size(), std::vector<int>(k_path_metrics.front().size()));
  std::vector<std::vector<int>> new_k_trellis_column(
      k_path_metrics.size(), std::vector<int>(k_path_metrics.front().size()));
  for (size_t i = 0; i < k_path_metrics.size(); i++) {
    // Identify which k-th path metric was used
    std::vector<std::vector<bool>> identifies(
        k_path_metrics.size(), std::vector<bool>(k_path_metrics.front().size(), false));
    // Loop L times
    for (size_t j = 0; j < k_path_metrics[i].size(); j++) {
      auto p = PathMetric(bits, k_path_metrics, identifies, i);
      new_k_path_metrics[i][j] = p.first;
      new_k_trellis_column[i][j] = p.second;
    }
  }

  k_path_metrics = new_k_path_metrics;
  trellis.push_back(new_k_trellis_column);
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

std::pair<int, int>
ViterbiCodec::PathMetric(const std::string &bits,
                         const std::vector<std::vector<int>> &prev_k_path_metrics,
                         std::vector<std::vector<bool>> &prev_identifies,
                         int state) const {
  // Find the source state that can transform to state
  int s = (state << 1) & ((1 << (constraint_ - 1)) - 1);
  int source_state1 = s | 0;
  int source_state2 = s | 1;

  // Find the first not using path metric
  size_t selected_index1 = 0;
  while (selected_index1 < prev_identifies[source_state1].size() && prev_identifies[source_state1][selected_index1]) {
    selected_index1++;
  }
  int k_path_metric1 = prev_k_path_metrics[source_state1][selected_index1];
  if (k_path_metric1 < std::numeric_limits<int>::max()) {
    k_path_metric1 += BranchMetric(bits, source_state1, state);
  }

  size_t selected_index2 = 0;
  while (selected_index2 < prev_identifies[source_state2].size() && prev_identifies[source_state2][selected_index2]) {
    selected_index2++;
  }
  int k_path_metric2 = prev_k_path_metrics[source_state2][selected_index2];
  if (k_path_metric2 < std::numeric_limits<int>::max()) {
    k_path_metric2 += BranchMetric(bits, source_state2, state);
  }

  k_path_metric1 <= k_path_metric2 ? prev_identifies[source_state1][selected_index1] = true : prev_identifies[source_state2][selected_index2] = true;

  return k_path_metric1 <= k_path_metric2 ? std::make_pair(k_path_metric1, source_state1) :
                                            std::make_pair(k_path_metric2, source_state2);
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
  for (size_t i = 0; i < x.size(); i++) {
    distance += x[i] != y[i];
  }
  return distance;
}

}// namespace lab