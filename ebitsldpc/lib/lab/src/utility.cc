#include "utility.h"

namespace lab {
namespace utility {
// using inline to avoid multi definition of function
void
MatrixProd(int *uu, int *cc, int **G, int dim, int len) {
  for (int n = 0; n < len; ++n) {
    cc[n] = 0;
  }
  for (int k = 0; k < dim; ++k) {
    for (int n = 0; n < len; ++n) {
      cc[n] ^= (uu[k] & G[k][n]);
    }
  }
}

void
ProbClip(double *xx, int len_xx) {
  for (int i = 0; i < len_xx; i++) {
    if (xx[i] < kSmallestProb) {
      xx[i] = kSmallestProb;
    } else if (xx[i] > 1.0 - kSmallestProb) {
      xx[i] = 1.0 - kSmallestProb;
    }
  }
}

std::ostream &
operator<<(std::ostream &os, const std::complex<double> &c) {
  os << std::fixed << std::setprecision(14)
     << std::resetiosflags(std::ios::showpos)
     << '('
     << std::showpos << c.real() << std::showpos << c.imag() << 'i'
     << ')'
     << std::resetiosflags(std::ios::showpos);
  return os;
}
}// namespace utility
}// namespace lab
