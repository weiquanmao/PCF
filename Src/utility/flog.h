#ifndef _LOG_UTILITY_H_FILE_
#define _LOG_UTILITY_H_FILE_

#include <iostream>
#include <fstream>

//#define USE_FLOG_NAMESPACE

#if defined(USE_FLOG_NAMESPACE)
namespace FLOG {
#endif

    bool isLogFileOpened();

    bool setDefaultLogFile();
    bool setLogFile(const char *logFile);
    const char* getLogFile();

    void setAppendMode();
    void setWriteMode();

    int flog(const char *format, ...);


    class FLogCout
    {
        friend bool _closeLog(bool cleanPath);
        friend bool _openLog(const char *logFile);
    public:
        FLogCout() { setDefaultLogFile(); }
        ~FLogCout() {
            if (_m_log_file_.is_open())
                _m_log_file_.close();
        }

        template<typename T>
        FLogCout& operator << (const T a)
        {
            std::cout << a;
            if (_m_log_file_.is_open()) {
                _m_log_file_ << a;
                _m_log_file_.flush();
            }
            return *this;
        }
    private:
        std::ofstream _m_log_file_;
    };

extern FLogCout clog;

#if defined(USE_FLOG_NAMESPACE)
}
#endif

#endif // !_LOG_UTILITY_H_FILE_
