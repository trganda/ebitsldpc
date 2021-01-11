#ifndef LAB_INCLUDE_VITERBI_H
#define LAB_INCLUDE_VITERBI_H

#include <vector>
#include <utility>
#include <string>
#include <cassert>
#include <limits>
#include <algorithm>

namespace lab {

class ViterbiCodec {
 public:
  // Note about Polynomial Descriptor of a Convolutional Encoder / Decoder.
  // A generator polymonial is built as follows: Build a binary number
  // representation by placing a 1 in each spot where a connection line from
  // the shift register feeds into the adder, and a zero elsewhere. There are 2
  // ways to arrange the bits:
  // 1. msb-current
  //    The MSB of the polynomial corresponds to the current input, while the
  //    LSB corresponds to the oldest input that still remains in the shift
  //    register.
  //    This representation is used by MATLAB. See
  //    http://radio.feld.cvut.cz/matlab/toolbox/comm/tutor124.html
  // 2. lsb-current
  //    The LSB of the polynomial corresponds to the current input, while the
  //    MSB corresponds to the oldest input that still remains in the shift
  //    register.
  //    This representation is used by the Spiral Viterbi Decoder Software
  //    Generator. See http://www.spiral.net/software/viterbi.html
  // We use 2.
  ViterbiCodec(int constraint, const std::vector<int> &polynomials);
  std::string Encode(const std::string &bits) const;
  std::string Decode(const std::string &bits) const;
  int constrain() const { return constraint_; }
  int num_parity_bits() const { return polynomials_.size(); }

 private:
  typedef std::vector<std::vector<int>> Trellis;
  typedef std::vector<std::vector<std::vector<int>>> PLV_Trellis;
  void initializeOutputs();
  // Given num_parity_bits() received bits, update path metrics of all states
  // in the current iteration, and append new traceback vector to trellis.
  void UpdatePathMetrics(const std::string &bits,
                         std::vector<int> &path_metrics,
                         Trellis &trellis) const;
  // Given num_parity_bits() received bits, compute and returns path
  // metric and its corresponding previous state.
  std::pair<int, int> PathMetric(const std::string &bits,
                                 const std::vector<int> &prev_path_metrics,
                                 int state) const;
  int BranchMetric(const std::string &bits,
                   int source_state,
                   int target_state) const;
  std::string Output(int current_state, int input) const;
 private:
  const int constraint_;
  const std::vector<int> polynomials_;

  // The output table.
  // The index is current input bit combined with previous inputs in the shift
  // register. The value is the output parity bits in string format for
  // convenience, e.g. "10". For example, suppose the shift register contains
  // 0b10 (= 2), and the current input is 0b1 (= 1), then the index is 0b110 (=
  // 6).
  std::vector<std::string> outputs_;
};

int ReverseBits(int num_bits, int input);

int HammingDistance(const std::string &x, const std::string &y);

}
#endif //LAB_INCLUDE_VITERBI_H
