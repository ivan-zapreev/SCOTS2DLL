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
                    const int m_ipt_dof_idx;
                    const int m_ss_size;
                    jdoubleArray m_args;
                protected:
                public:

                    /**
                     * Constructs a new controller wrapper for the java class
                     * @param env to jni environment
                     * @param person_class the class to construct the wrapper for
                     * @param name the individual class name
                     * @param ipt_dof_idx the input dof index
                     * @param ss_size the input state space size
                     */
                    ctrl_wrapper(JNIEnv * const env, const jclass & person_class,
                            const string name, const jint ipt_dof_idx, const int ss_size)
                    : m_env(env), m_person_class(person_class),
                    m_eval_mthd(m_env->GetStaticMethodID(m_person_class, "evaluate_0", "([D)D")),
                    m_name(name), m_ipt_dof_idx(ipt_dof_idx), m_ss_size(ss_size),
                    m_args(m_env->NewDoubleArray(m_ss_size)) {
                    }

                    /**
                     * Allows to get the input dof index
                     * @return the input dof index
                     */
                    inline const int & get_dof_idx() const {
                        return m_ipt_dof_idx;
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
                    inline double compute_input(const double *state) {
                        m_env->SetDoubleArrayRegion(m_args, 0, m_ss_size, state);
                        return m_env->CallStaticDoubleMethod(
                                m_person_class, m_eval_mthd, m_args);
                    }
                };
            }
        }
    }
}

#endif /* CTRL_WRAPPER_HH */

