#include "mat.h"

namespace lab {
#ifdef USE_MATIO

Mat::Mat(std::string &filename)
    : filename_(filename), matfp_(nullptr) {}

Mat::~Mat() {
  if (nullptr != matfp_) {
    Mat_Close(matfp_);
  }
}

void
Mat::Open() {
  matfp_ = Mat_CreateVer(
      filename_.c_str(), nullptr,
      MAT_FT_DEFAULT);
  if (nullptr == matfp_) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Creating file failed, file name is "
                             << filename_ << std::endl;
    exit(-1);
  }
}

void
Mat::Close() {
  if (nullptr != matfp_) {
    Mat_Close(matfp_);
    matfp_ = nullptr;
  }
}

void
Mat::WriteInt(
    const std::string &varname,
    int32_t data) {
  if (nullptr == matfp_) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Create of Open the file first" << std::endl;
    exit(-1);
  }
  size_t dims[2] = {1, 1};
  int32_t datas[1] = {data};
  matvar_t *mat_var = Mat_VarCreate(
      varname.c_str(),
      MAT_C_INT32, MAT_T_INT32,
      2, dims, datas, 0);
  if (nullptr == matfp_) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Error creating variable for " << varname << std::endl;
    exit(-1);
  } else {
    Mat_VarWrite(matfp_, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);
  }
  logger::Log::get().set_log_level(logger::Info);
  LOG(logger::Info, false) << "Writed " << varname << " to " << filename_ << std::endl;
}

void
Mat::WriteDouble(const std::string &varname, double data) {
  if (nullptr == matfp_) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Create of Open the file first" << std::endl;
    exit(-1);
  }
  size_t dims[2] = {1, 1};
  double datas[1] = {data};
  matvar_t *mat_var = Mat_VarCreate(
      varname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, datas, 0);
  if (nullptr == mat_var) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Error creating variable for " << varname << std::endl;
    exit(-1);
  } else {
    Mat_VarWrite(matfp_, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);
  }
  logger::Log::get().set_log_level(logger::Info);
  LOG(logger::Info, false) << "Writed " << varname << " to " << filename_ << std::endl;
}

void
Mat::WriteComplex(const std::string &varname, std::complex<double> data) {
  if (nullptr == matfp_) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Create of Open the file first" << std::endl;
    exit(-1);
  }
  size_t dims[2] = {1, 1};
  double x[1] = {data.real()},
      y[1] = {data.imag()};
  struct mat_complex_split_t datas = {x, y};
  matvar_t *mat_var = Mat_VarCreate(
      varname.c_str(),
      MAT_C_DOUBLE, MAT_T_DOUBLE,
      2, dims, &datas, MAT_F_COMPLEX);
  if (nullptr == mat_var) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Error creating variable for " << varname << std::endl;
    exit(-1);
  } else {
    Mat_VarWrite(matfp_, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);
  }
  logger::Log::get().set_log_level(logger::Info);
  LOG(logger::Info, false) << "Writed " << varname << " to " << filename_ << std::endl;
}

void
Mat::WriteVector(const std::string &varname, const std::vector<int32_t> &data) {
  if (nullptr == matfp_) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Create of Open the file first" << std::endl;
    exit(-1);
  }
  size_t dims[2] = {data.size(), 1};
  int32_t *data_wrap = new int32_t[data.size()];
  for (size_t i = 0; i < data.size(); i++) {
    data_wrap[i] = data[i];
  }
  matvar_t *mat_var = Mat_VarCreate(varname.c_str(), MAT_C_INT32, MAT_T_INT32, 2, dims, data_wrap, 0);
  delete[] data_wrap;
  if (nullptr == mat_var) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Error creating variable for " << varname << std::endl;
    exit(-1);
  } else {
    Mat_VarWrite(matfp_, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);
  }
  logger::Log::get().set_log_level(logger::Info);
  LOG(logger::Info, false) << "Writed " << varname << " to " << filename_ << std::endl;
}

void
Mat::WriteVector(const std::string &varname, const std::vector<double> &data) {
  if (nullptr == matfp_) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Create of Open the file first" << std::endl;
    exit(-1);
  }
  size_t dims[2] = {data.size(), 1};
  int32_t *data_wrap = new int32_t[data.size()];
  for (size_t i = 0; i < data.size(); i++) {
    data_wrap[i] = data[i];
  }
  matvar_t *mat_var = Mat_VarCreate(
      varname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, data_wrap, 0);
  delete[] data_wrap;
  if (nullptr == mat_var) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Error creating variable for " << varname << std::endl;
    exit(-1);
  } else {
    Mat_VarWrite(matfp_, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);
  }
  logger::Log::get().set_log_level(logger::Info);
  LOG(logger::Info, false) << "Writed " << varname << " to " << filename_ << std::endl;
}

void
Mat::WriteVector(const std::string &varname, const std::vector<std::complex<double>> &data) {
  if (nullptr == matfp_) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Create of Open the file first" << std::endl;
    exit(-1);
  }
  size_t dims[2] = {data.size(), 1};
  double x[data.size()],
      y[data.size()];
  for (size_t i = 0; i < data.size(); i++) {
    x[i] = data[i].real();
    y[i] = data[i].imag();
  }
  struct mat_complex_split_t datas = {x, y};
  matvar_t *mat_var = Mat_VarCreate(
      varname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &datas,
      MAT_F_COMPLEX);
  if (nullptr == mat_var) {
    logger::Log::get().set_log_level(logger::Error);
    LOG(logger::Error, true) << "Error creating variable for " << varname << std::endl;
    exit(-1);
  } else {
    Mat_VarWrite(matfp_, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);
  }
  logger::Log::get().set_log_level(logger::Info);
  LOG(logger::Info, false) << "Writed " << varname << " to " << filename_ << std::endl;
}

#endif
}// namespace lab
