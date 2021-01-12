#ifndef LAB_INCLUDE_VITERBI_H
#define LAB_INCLUDE_VITERBI_H

#include <algorithm>
#include <cassert>
#include <limits>
#include <string>
#include <utility>
#include <vector>

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
  std::vector<std::string> PLVDecode(const std::string &bits, unsigned int L) const;
  int constrain() const { return constraint_; }
  int num_parity_bits() const { return polynomials_.size(); }

 private:
  // The trellis table.
  // The trellis[i][j] is mean that the prev state occupied by the best path to
  // state j at instant i. for a constraint 3, and polynomial [7, 5] coder, decode
  // "0011100001100111111000101100111011", the transpose of trellis is follow,
  //  instant  0  ******************************************* 16
  //  state 0  0  0  0  1  0  1  1  1  0  0  1  0  1  0  0  0  1
  //  state 1  2  2  2  3  2  3  3  3  2  2  3  2  3  2  2  2  3
  //  state 2  0  0  0  1  0  1  1  1  0  0  1  0  1  0  0  0  1
  //  state 3  2  2  2  3  2  3  3  3  2  2  3  2  3  2  2  2  3
  // the last state occupied by best path is 0 that calculate from path_metric,
  // Then we can get the best path follow,
  //  instant  0  *************************************************************************** 16
  //  state 0  0 -  0    0    1    0    1    1 -> 1 -  0    0    1    0 -> 1 -> 0 -  0    0 -> 1
  //  state 1  2 -  2 -> 2 -  3    2    3 -> 3 -  3 -  2 -> 2 -  3 -> 2 -  3    2 -  2 -> 2 -  3
  //  state 2  0 -> 0 -  0 -> 1 -  0    1 -  1    1 -> 0 -  0 -> 1 -  0    1    0 -> 0 -  0    1
  //  state 3  2    2    2    3 -> 2 -> 3 -  3    3    2    2    3    2    3    2    2    2    3
  // the states occupied by best path was "02123310212100210", and we can get
  // the decoded information (removed flushing bits) is "01011100101000100"
  typedef std::vector<std::vector<int>> Trellis;
  // The Parallel List Viterbi trellis table
  // 3 dimension:[i][j][k]
  //   1st-dim: time-instant i
  //   2ed-dim: state j
  //   3rd-dim: traceback for k-best path
  typedef std::vector<std::vector<std::vector<int>>> PLV_Trellis;
  void initializeOutputs();
  // Given num_parity_bits() received bits, update path metrics of all states
  // in the current iteration, and append new traceback vector to trellis.
  void UpdatePathMetrics(const std::string &bits,
                         std::vector<int> &path_metrics,
                         Trellis &trellis) const;
  void UpdatePathMetrics(const std::string &bits,
                         std::vector<std::vector<int>> &k_path_metrics,
                         PLV_Trellis &trellis) const;
  // Given num_parity_bits() received bits, compute and returns path
  // metric and its corresponding previous state.
  std::pair<int, int> PathMetric(const std::string &bits,
                                 const std::vector<int> &prev_path_metrics,
                                 int state) const;
  std::pair<int, int> PathMetric(const std::string &bits,
                                 const std::vector<std::vector<int>> &prev_k_path_metrics,
                                 std::vector<std::vector<bool>> &prev_identifies,
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

int
ReverseBits(int num_bits, int input);

int
HammingDistance(const std::string &x, const std::string &y);

}// namespace lab
#endif//LAB_INCLUDE_VITERBI_H
