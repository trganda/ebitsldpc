#ifndef LAB_LOG_H
#define LAB_LOG_H

#include <chrono>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>

#ifndef __FILENAME__
#define __FILENAME__ __FILE__
#endif
#define LOG(level, flag)                                                                 \
  if (level <= lab::logger::Log::get().log_level())                                      \
  lab::logger::Log::get().log_stream(flag) << '[' << lab::logger::Log::get_time() << ']' \
                                           << '[' << __FILENAME__ << ':' << std::dec << __LINE__ << "] "
namespace lab::logger {
enum Level {
  Error,
  Info,
};
const std::map<Level, std::string> colored{
    {Error, " \x1b[31;1m[ERROR]\x1b[0m "},
    {Info, " \x1b[32;1m[INFO]\x1b[0m "}};
//Courtesy of http://wordaligned.org/articles/cpp-streambufs#toctee-streams
class TeeBuf : public std::streambuf {
 public:
  // Construct a streambuf which tees output to both input
  // streambufs.
  explicit TeeBuf(std::streambuf *sbtofile, std::streambuf *sbtoscreen);
  void set_flag(bool flag);

 private:
  // This tee buffer has no buffer. So every character "overflows"
  // and can be put directly into the teed buffers.
  int overflow(int c) override;
  // Sync both teed buffers.
  int sync() override;

 private:
  std::streambuf *strbuf_to_file_;
  std::streambuf *strbuf_to_stdout_;
  bool flag_;
};
class TeeStream : public std::ostream {
 public:
  explicit TeeStream(std::ostream &o1, std::ostream &o2);
  TeeBuf &tbuf();

 private:
  TeeBuf tbuf_;
};
class Log {
 public:
  ~Log() = default;
  void set_log_stream(std::ostream &stream);
  Log &set_log_level(Level level);
  Level log_level() const;
  std::ostream &log_stream();
  std::ostream &log_stream(bool flag);
  static Log &get();
  static std::string get_time();
  void log(const std::string &message, Level level, bool flag);

 private:
  void log(const std::string &message, bool flag);

 private:
  Level log_level_;
  std::ostream *log_stream_;
  std::mutex log_mutex;
};

inline void
ERROR(const std::string &message, bool both_to_stdout) {
  Log::get().log(message, Level::Error, both_to_stdout);
}

inline void
INFO(const std::string &message, bool both_to_stdout) {
  Log::get().log(message, Level::Info, both_to_stdout);
}
}// namespace lab::logger
#endif
