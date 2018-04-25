/*
 * File:   scots2dll.cc
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
 * Created on January 14, 2018, 8:36 PM
 */

#ifndef JNI_THROW_HH
#define JNI_THROW_HH

#include <jni.h>
#include <string>

#include "info_logger.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                const char * const FileNotFoundException = "java/io/FileNotFoundException";
                const char * const IllegalStateException = "java/lang/IllegalStateException";
                const char * const IllegalArgumentException = "java/lang/IllegalArgumentException";
                const char * const ClassNotFoundException = "java/lang/ClassNotFoundException";

                jint throwException(info_logger &logger, JNIEnv *env,
                        const char * const class_name, const char * const message) {
                    LOG(logger, "ERROR (" << class_name << "): " << message);
                    jclass ex_class = env->FindClass(class_name);
                    if (ex_class == NULL) {
                        return -1;
                    } else {
                        return env->ThrowNew(ex_class, message);
                    }
                }

            }
        }
    }
}

#endif /* JNI_THROW_HH */

