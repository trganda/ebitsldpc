#include <chrono>
#include <iostream>

#include "log.h"
#include "randnum.h"
#include "toml.hpp"
#include "simulator.h"

int
main() {
  // Start measuring time
  auto begin = std::chrono::high_resolution_clock::now();

  std::string dirname = "logs";
  std::string logFileName = lab::logger::Log::get_time() + "-kmldpc.logger";
  std::ofstream logFile(dirname + "/" + logFileName);
  lab::logger::TeeStream logTee(logFile, std::cout);
  if (logFile.is_open() && logFile.good())
    lab::logger::Log::get().set_log_stream(logTee);
  else
    lab::logger::Log::get().set_log_stream(std::cout);
  int flag0, flag1;
  flag0 = 0;
  lab::CLCRandNum::Get().SetSeed(flag0);
  flag1 = 0;
  lab::CWHRandNum::Get().SetSeed(flag1);
  lab::logger::Log::get().set_log_level(lab::logger::Info);
  lab::logger::INFO("Start simulation", true);
  std::ifstream ifs("config.toml", std::ios_base::binary);
  if (ifs.is_open()) {
    auto arguments = toml::parse(ifs);
    int constraint = 3;
    std::vector<int> polynomials = {7, 5};
    lab::ViterbiCodec ccode (constraint, polynomials);
    auto ret = ccode.Decode("0011100001100111111000101100111011");
    ifs.close();
  } else {
    lab::logger::ERROR("Encouter error while opening config.toml", true);
  }
  lab::logger::INFO("Simulation done", true);

  // Stop measuring time and calculate the elapsed time
  auto end = std::chrono::high_resolution_clock::now();
  auto total_time_minutes = std::chrono::duration_cast<std::chrono::minutes>(end - begin).count();
  auto total_time_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
  auto total_time_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

  total_time_milliseconds -= total_time_minutes * 60 * 1000;
  total_time_seconds -= total_time_minutes * 60;

  std::stringstream stream;
  stream << "Total time cost: "
         << total_time_minutes << "min:"
         << total_time_seconds << "sec:"
         << total_time_milliseconds << "ms";
  lab::logger::INFO(stream.str(), true);
  return 0;
}