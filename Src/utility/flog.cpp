#include "flog.h"
#include <stdarg.h>
#include <stdio.h>
#if defined(USE_FLOG_NAMESPACE)
namespace FLOG {
#endif
    FLogCout clog;
    FILE *_log_file_ = 0;
    char _log_file_path_[1024];
    bool _b_log_file_opened_ = false;
    char _log_mode_ = 'w';
    const char _log_file_path_default_[] = "../flog.txt";
    
    bool _closeLog(bool cleanPath)
    {
        bool bRet = true;
        if (_b_log_file_opened_) {
            if (fclose(_log_file_) == 0) {
                _log_file_ = 0;
                clog._m_log_file_ = std::ofstream(_log_file_);
                _b_log_file_opened_ = false;
                if (cleanPath)
                     sprintf(_log_file_path_, "");
            }
            else {
                printf("[FLOG ERROR]: Failed to close log file %s.", _log_file_path_);
                bRet = false;
            }
        }
        return bRet;
    }
    bool _openLog(const char *logFile)
    {
        bool bRet = true;
        if (logFile == 0)
            bRet = _closeLog(true);
        else {
            if (_b_log_file_opened_ && !_closeLog(false))
                return false;

            FILE *log_file_temp = fopen(logFile, &_log_mode_);
            if (log_file_temp != 0) {               
                _log_file_ = log_file_temp;
                clog._m_log_file_ = std::ofstream(_log_file_);
                sprintf(_log_file_path_, logFile);
                _b_log_file_opened_ = true;
                bRet = true;
            }
            else {
                printf("[FLOG ERROR]: Failed to open log file %s.", _log_file_path_);
                bRet = false;
            }
        }

        return bRet;
    }
    
    bool isLogFileOpened() { return _b_log_file_opened_; }
    bool setDefaultLogFile() { return _openLog(_log_file_path_default_); }
    bool setLogFile(const char *logFile) { return _openLog(logFile); }
    const char* getLogFile() { return _log_file_path_; }

    void setAppendMode() { _log_mode_ = 'a'; }
    void setWriteMode() { _log_mode_ = 'w'; }

    int flog(const char *format, ...)
    {        
        va_list ap1;
        va_start(ap1, format);
        int n1 = vprintf(format, ap1);
        va_end(ap1);

        if (!_b_log_file_opened_)
            _openLog(_log_file_path_default_);

        if (_b_log_file_opened_) {
            va_list ap2;
            va_start(ap2, format);
            vfprintf(_log_file_, format, ap2);
            va_end(ap2);
            fflush(_log_file_);
        }
        return n1;
    }

#if defined(USE_FLOG_NAMESPACE)
}
#endif