/* 
 * File:   info_logger.hh
 * Author: Dr. Ivan S. Zapreev
 *
 * Visit my Linked-in profile:
 *      <https://nl.linkedin.com/in/zapreevis>
 * Visit my GitHub:
 *      <https://github.com/ivan-zapreev>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.#
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Created on April 23, 2018, 10:50 AM
 */

#ifndef INFO_LOGGER_HH
#define INFO_LOGGER_HH

#include <mutex>
#include <iostream>

using namespace std;

#define LOG(MSG) \
                { \
                    time_t now = time(0);\
                    string time(ctime(&now));\
                    time = time.substr(0, time.length() - 1);\
                    {\
                        info_logger & log = get_log();\
                        lock_guard<recursive_mutex> guard(log.m_mutex); \
                        log.m_log_file << time << " " << __FILE__ << ":" \
                        << __LINE__ << " > " << MSG << std::flush << std::endl; \
                    }\
                }\

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * This class is used for basic logging from the native code, it is not thread safe
                 */
                class info_logger {
                public:
                    //Stores the recursive mutex for thread safety
                    recursive_mutex m_mutex;

                    //Stores the log file
                    ofstream m_log_file;

                    info_logger() : m_log_file() {
                    }

                    virtual ~info_logger() {
                        close_log();
                    }

                    void open_log(const char * file_name) {
                        //Read controller from file
                        string log_file_name = string(file_name) + ".sr.log";

                        //Open the file in a thread-safe way
                        {
                            lock_guard<recursive_mutex> guard(m_mutex);

                            //Stop logging if the logging was started
                            close_log();

                            //Open the file for writing and clear it
                            m_log_file.open(log_file_name);
                            m_log_file.clear();
                        }

                        m_log_file << "Opened the log-file: " << log_file_name;
                    }

                    void close_log() {
                        //Close the file in the thread-safe manner
                        lock_guard<recursive_mutex> guard(m_mutex);

                        //Close the file if it is open
                        if (m_log_file.is_open()) {
                            m_log_file << std::flush;
                            m_log_file.close();
                        }
                    }
                };
                
                //Stores the info logger object
                static info_logger m_log;

                /**
                 * Allows to (re-)open the info log
                 * @param file_name the log file name
                 */
                inline void open_info_log(const char * file_name) {
                    //Open logging
                    m_log.open_log(file_name);
                }

                /**
                 * Allows to get the info log object
                 * @return the info log object
                 */
                inline info_logger & get_log() {
                    return m_log;
                }

                /**
                 * Allows to close the info log, if it is open
                 */
                inline void close_info_log() {
                    m_log.close_log();
                }
                
            }
        }
    }
}

#endif /* INFO_LOGGER_HH */

