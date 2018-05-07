/*
 * File:   ctrl_wrapper.hh
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
 * Created on January 14, 2018, 6:49 PM
 */

#ifndef CTRL_WRAPPER_HH
#define CTRL_WRAPPER_HH

#include <jni.h>
#include <string>

#include "jni_throw.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * This class represents a generic CTRL wrapper
                 */
                class ctrl_wrapper {
                private:
                    JNIEnv * const m_env;
                    const jclass m_person_class;
                    const jmethodID m_eval_mthd;
                    const string m_name;
                    const int m_ss_dim;
                    const int m_is_dim;
                    jdoubleArray m_state;
                    jdoubleArray m_input;
                protected:
                public:

                    /**
                     * Constructs a new controller wrapper for the java class
                     * @param env to jni environment
                     * @param person_class the class to construct the wrapper for
                     * @param name the individual class name
                     * @param ss_size the number of state space dimensions
                     * @param is_size the number of input space dimensions
                     */
                    ctrl_wrapper(JNIEnv * const env, const jclass & person_class,
                            const string name, const int ss_dim, const int is_dim)
                    : m_env(env), m_person_class(person_class),
                    m_eval_mthd(env->GetStaticMethodID(m_person_class, "evaluate", "([D[D)V")),
                    m_name(name), m_ss_dim(ss_dim), m_is_dim(is_dim),
                    m_state(env->NewDoubleArray(m_ss_dim)),
                    m_input(env->NewDoubleArray(m_is_dim)){
                    }
                    
                    virtual ~ctrl_wrapper() {
                    }
                    
                    /**
                     * Allows to get the state-space dimensionality
                     * @return the state-space dimensionality
                     */
                    inline int get_ss_dim() {
                        return m_ss_dim;
                    }
                    
                    /**
                     * Allows to get the input-space dimensionality
                     * @return the input-space dimensionality
                     */
                    inline int get_is_dim() {
                        return m_is_dim;
                    }

                    /**
                     * Allows to get the wrapper's name
                     * @return the wrapper's name
                     */
                    inline const string & get_name() const {
                        return m_name;
                    }

                    /**
                     * Computes the input value for the given state
                     * @param state the state to call the method for
                     * @return the computed input value for this state
                     */
                    inline void compute_input(const double *state, double * input) {
                        m_env->SetDoubleArrayRegion(m_state, 0, m_ss_dim, state);
                        m_env->CallVoidMethod(m_person_class, m_eval_mthd, m_state, m_input);
                        m_env->GetDoubleArrayRegion(m_input, 0, m_is_dim, input);
                    }
                };
            }
        }
    }
}

#endif /* CTRL_WRAPPER_HH */

