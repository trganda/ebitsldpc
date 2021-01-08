#ifndef LAB_RAND_NUM_H
#define LAB_RAND_NUM_H

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

namespace lab {
class RandNum {
 public:
  explicit RandNum() = default;
  virtual ~RandNum() = default;
  virtual void SetSeed(int flag) = 0;
  virtual void PrintState(FILE *fp) = 0;
  virtual double Uniform() = 0;
  virtual void Normal(double *, int) = 0;
};
/*
 * The following generator employs the linear-congruential method,
 * and specifically uses a choice of multiplier that was proposed
 * as a standard by Stephen K. Park et al. in "Technical correspondence,"
 * Communications of the ACM36(1993), number 7, 108-110
 */

class CLCRandNum : public RandNum {
 public:
  explicit CLCRandNum();
  ~CLCRandNum() override = default;
  void SetSeed(int flag) final;
  void PrintState(FILE *fp) final;
  double Uniform() final;
  void Normal(double *nn, int len_nn) final;
  void Normal(std::complex<double> &nn);
  void Normal(std::vector<std::complex<double>> &nn);
  // Singleton pattern
  static CLCRandNum &Get();

 private:
  long int state;
  int A;
  long M;
  int Q;
  int R;
};
/*
 *The following generator employs the Wichman-Hill algorithm
 */
class CWHRandNum : public RandNum {
 public:
  explicit CWHRandNum();
  ~CWHRandNum() override = default;
  void SetSeed(int flag) final;
  void PrintState(FILE *fp) final;
  double Uniform() final;
  void Normal(double *nn, int len_nn) final;
  static CWHRandNum &Get();

 private:
  int X;
  int Y;
  int Z;
};
}// namespace lab
#endif// RAND_NUM_H
