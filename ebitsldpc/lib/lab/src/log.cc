#include "log.h"

namespace lab::logger {
TeeBuf::TeeBuf(std::streambuf *sbtofile, std::streambuf *sbtoscreen)
    : strbuf_to_file_(sbtofile),
      strbuf_to_stdout_(sbtoscreen),
      flag_(true) {}

void
TeeBuf::set_flag(bool flag) {
  this->flag_ = flag;
}

int
TeeBuf::overflow(int c) {
  if (c == EOF) {
    return !EOF;
  } else {
    const int r1 = strbuf_to_file_->sputc(c);
    if (flag_) {
      const int r2 = strbuf_to_stdout_->sputc(c);
      return r1 == EOF || r2 == EOF ? EOF : c;
    } else {
      return r1 == EOF ? EOF : c;
    }
  }
}

int
TeeBuf::sync() {
  const int r1 = strbuf_to_file_->pubsync();
  if (flag_) {
    const int r2 = strbuf_to_stdout_->pubsync();
    return r1 == 0 && r2 == 0 ? 0 : -1;
  } else {
    return r1 == 0 ? 0 : -1;
  }
}

TeeStream::TeeStream(std::ostream &o1, std::ostream &o2)
    : std::ostream(&tbuf_),
      tbuf_(o1.rdbuf(), o2.rdbuf()) {}

TeeBuf &
TeeStream::tbuf() {
  return tbuf_;
}

void
Log::set_log_stream(std::ostream &stream) {
  log_stream_ = &stream;
}

Log &
Log::set_log_level(Level level) {
  log_level_ = level;
  return *this;
}

Level
Log::log_level() const {
  return log_level_;
}

std::ostream &
Log::log_stream() {
  return *log_stream_;
}

std::ostream &
Log::log_stream(bool flag) {
  auto ptr = (TeeStream *) log_stream_;
  ptr->tbuf().set_flag(flag);
  return *log_stream_;
}

std::string
Log::get_time() {
  std::time_t secSinceEpoch = std::chrono::system_clock::to_time_t(
      std::chrono::system_clock::now());
  struct tm *calendarTime = localtime(&secSinceEpoch);
  char usrdefFormat[50] = {0};
  strftime(
      usrdefFormat, 50,
      "%Y-%m-%d %H:%M:%S", calendarTime);
  return std::string(usrdefFormat);
}

void
Log::log(const std::string &message, Level level, bool flag) {
  if (level > this->log_level_) {
    return;
  }
  std::stringstream stream;
  stream << '[' << get_time() << ']'
         << colored.at(level)
         << message << std::endl;
  log(stream.str(), flag);
}

void
Log::log(const std::string &message, bool flag) {
  std::lock_guard<std::mutex> lock(this->log_mutex);
  log_stream(flag);
  log_stream() << message;
}

Log &
Log::get() {
  static Log instance;
  return instance;
}
}// namespace lab::logger
