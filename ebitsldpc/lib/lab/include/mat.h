#ifndef LAB_MAT_H
#define LAB_MAT_H

#include <complex>
#include <string>
#include <vector>

#ifdef USE_MATIO
#include <matio.h>
#endif

#include "log.h"

namespace lab {

// Encapsulation of matio https://github.com/tbeu/matio
#ifdef USE_MATIO
class Mat {
 public:
    explicit Mat(std::string &filename);
    ~Mat();

 public:
    void Open();
    void Close();
    void WriteInt(const std::string &varname, int32_t data);
    void WriteDouble(const std::string &varname, double data);
    void WriteComplex(const std::string &varname, std::complex<double> data);
    void WriteVector(const std::string &varname, const std::vector<int32_t> &data);
    void WriteVector(const std::string &varname, const std::vector<double> &data);
    void WriteVector(const std::string &varname, const std::vector<std::complex<double>> &data);

 private:
    std::string filename_;
    mat_t *matfp_;
};
#endif
}// namespace lab
#endif
