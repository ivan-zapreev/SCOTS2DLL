/* 
 * File:   config.hh
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
 * Created on April 23, 2018, 11:26 AM
 */

#ifndef CONFIG_WRAPPER_HH
#define CONFIG_WRAPPER_HH

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * The fitness types enumeration
                 */
                enum fitness_type {
                    UNDEF_FITNESS = 0,
                    EXACT_FITNESS = 1,
                    ATANG_FITNESS = 2,
                    INVER_FITNESS = 3
                };

                //Stores the maximum attractor size
                constexpr double MAX_ATTRACTOR_SIZE = 0.5;
                constexpr double EPSILON = 1e-2;
                //Stores the maximum stare deviation
                constexpr double MAX_STATE_DEVIATION = (MAX_ATTRACTOR_SIZE - EPSILON);

                /**
                 * Represents the configuration object class swapper
                 */
                class config_wrapper {
                public:

                    //The number of state-space dimensions
                    const int m_num_ss_dim;
                    //The number of input-space dimensions
                    const int m_num_is_dim;
                    //Stores the attractor size
                    const double m_attr_size;
                    //Stores the scaling flag
                    const bool m_is_scale;
                    //Stores the fitness type
                    const fitness_type m_ftn_type;
                    //Stores the flag for complex fitness
                    const bool m_is_complex;
                    //Stores the fitness function scaling factor for the complex fitness types INVERSE and ATANGENT
                    const double m_ftn_scale;
                    //True if the Monte-Carlo Fitness computations are to be used
                    const bool m_is_mc;
                    //True if the Recursive Stratified Sampling (RSS) is to be used
                    const bool m_is_rss;
                    //Stores the initial sample size for the Monte Carlo sampling
                    const long m_sample_size;
                    //Stores the minimum sample size the RSS is applied to
                    const long m_bisect_size;
                    //Stores the fraction of sample to be used for bisection
                    const double m_bisect_ratio;

                    /**
                     * The default constructor
                     */
                    config_wrapper()
                    : m_num_ss_dim(0), m_num_is_dim(0), m_attr_size(0.0),
                    m_is_scale(false), m_ftn_type(fitness_type::EXACT_FITNESS),
                    m_is_complex(false), m_ftn_scale(1.0), m_is_mc(false),
                    m_is_rss(false), m_sample_size(0), m_bisect_size(0),
                    m_bisect_ratio(0.0) {
                    }

                    /**
                     * The basic copy constructor
                     * @param other the object to copy from
                     */
                    config_wrapper(const config_wrapper & other)
                    : m_num_ss_dim(other.m_num_ss_dim), m_num_is_dim(other.m_num_is_dim),
                    m_attr_size(other.m_attr_size), m_is_scale(other.m_is_scale),
                    m_ftn_type(other.m_ftn_type), m_is_complex(other.m_is_complex),
                    m_ftn_scale(other.m_ftn_scale), m_is_mc(other.m_is_mc),
                    m_is_rss(other.m_is_rss), m_sample_size(other.m_sample_size),
                    m_bisect_size(other.m_bisect_size),
                    m_bisect_ratio(other.m_bisect_ratio) {
                    }

                    /**
                     * The basic copy assignment operator
                     * @param other the object to copy from
                     * @return the changed object assigned to
                     */
                    config_wrapper & operator=(const config_wrapper & other) {
                        const_cast<int&> (m_num_ss_dim) = other.m_num_ss_dim;
                        const_cast<int&> (m_num_is_dim) = other.m_num_is_dim;
                        const_cast<double&> (m_attr_size) = other.m_attr_size;
                        const_cast<bool&> (m_is_scale) = other.m_is_scale;
                        const_cast<fitness_type&> (m_ftn_type) = other.m_ftn_type;
                        const_cast<bool&> (m_is_complex) = other.m_is_complex;
                        const_cast<double&> (m_ftn_scale) = other.m_ftn_scale;
                        const_cast<bool&> (m_is_mc) = other.m_is_mc;
                        const_cast<bool&> (m_is_rss) = other.m_is_rss;
                        const_cast<long&> (m_sample_size) = other.m_sample_size;
                        const_cast<long&> (m_bisect_size) = other.m_bisect_size;
                        const_cast<double&> (m_bisect_ratio) = other.m_bisect_ratio;

                        return *this;
                    }

                    /**
                     * The basic constructor
                     * @param env the JNI environment
                     * @param obj the configuration object to wrap
                     * @param total_dim the total number of controller dimensions
                     */
                    config_wrapper(JNIEnv * env, jobject obj, const int total_dim)
                    : m_num_ss_dim(get_int_field(env, obj, "m_ss_size")),
                    m_num_is_dim(total_dim - m_num_ss_dim),
                    m_ftn_type((fitness_type) get_int_field(env, obj, "m_ftn_type")),
                    m_attr_size(get_double_field(env, obj, "m_attr_size")),
                    m_is_scale(get_bool_field(env, obj, "m_is_scale")),
                    m_is_complex(get_bool_field(env, obj, "m_is_complex")),
                    m_ftn_scale(get_ftn_scale(env, obj)),
                    m_is_mc(get_bool_field(env, obj, "m_is_mc")),
                    m_is_rss(get_bool_field(env, obj, "m_is_rss")),
                    m_sample_size(get_long_field(env, obj, "m_sample_size")),
                    m_bisect_size(get_long_field(env, obj, "m_bisect_size")),
                    m_bisect_ratio(get_double_field(env, obj, "m_bisect_ratio")) {
                    }

                    /**
                     * Allows to verify the consistency of the data, if the 
                     * data is not consistent then an exception is thrown
                     * @patam log the info logger class
                     * @param env the JNI environment
                     */
                    void finalize(info_logger & log, JNIEnv * env) const {
                        //Check the state space size
                        if (m_num_ss_dim < 1 || m_num_is_dim < 1) {
                            (void) throwException(log, env, IllegalArgumentException,
                                    "The state-space size must be [1, #dofs)");
                        }

                        //Check the attractor value
                        if (m_is_complex) {
                            if (m_attr_size < 0 || m_attr_size >= MAX_ATTRACTOR_SIZE) {
                                (void) throwException(log, env, IllegalArgumentException,
                                        "Improper attractor value, must be within [0, 0.5)!");
                            }
                        }

                        //If we have Monte Carlo fitness
                        if (m_is_mc) {
                            //Check on the sample size
                            if (m_sample_size < 1) {
                                (void) throwException(log, env, IllegalArgumentException,
                                        "The sample size for Monte Carlo simulations must be >= 1!");
                            }
                            //If we have Recursive Stratified Sampling
                            if (m_is_rss) {
                                //Check on the bisection size
                                if (m_bisect_size < 1) {
                                    (void) throwException(log, env, IllegalArgumentException,
                                            "The minimum bisection sample size for recursive stratified sampling must be >= 1!");
                                }
                                //Check on the bisection ratio
                                if ((m_bisect_ratio <= 0.0) || (m_bisect_ratio >= 1.0)) {
                                    (void) throwException(log, env, IllegalArgumentException,
                                            "The bisection ratio for recursive stratified sampling must be within (0.0,1.0)!");
                                }
                            }
                        }
                    }

                private:

                    /**
                     * Allows to get fitness scale, if the configuration 
                     * value is 0.0 then it is re-set to 1.0
                     * @param env the JNI environment
                     * @param obj the object itself
                     * @return the requested fitness scale value
                     */
                    double get_ftn_scale(JNIEnv * env, jobject obj) const {
                        const double ftn_scale = get_double_field(env, obj, "m_ftn_scale");
                        return (ftn_scale == 0.0) ? 1.0 : ftn_scale;
                    }

                    /**
                     * Allows to get the Java object field value
                     * @param env the JNI environment
                     * @param obj the object itself
                     * @param fid_name the field name
                     * @return the field value
                     */
                    int get_int_field(JNIEnv * env, jobject obj,
                            const char * fid_name) const {
                        jclass cls = env->GetObjectClass(obj);
                        jfieldID fid = env->GetFieldID(cls, fid_name, "I");
                        return env->GetIntField(obj, fid);
                    }

                    /**
                     * Allows to get the Java object field value
                     * @param env the JNI environment
                     * @param obj the object itself
                     * @param fid_name the field name
                     * @return the field value
                     */
                    long get_long_field(JNIEnv * env, jobject obj,
                            const char * fid_name) const {
                        jclass cls = env->GetObjectClass(obj);
                        jfieldID fid = env->GetFieldID(cls, fid_name, "J");
                        return env->GetIntField(obj, fid);
                    }

                    /**
                     * Allows to get the Java object field value
                     * @param env the JNI environment
                     * @param obj the object itself
                     * @param fid_name the field name
                     * @return the field value
                     */
                    double get_double_field(JNIEnv * env, jobject obj,
                            const char * fid_name) const {
                        jclass cls = env->GetObjectClass(obj);
                        jfieldID fid = env->GetFieldID(cls, fid_name, "D");
                        return env->GetDoubleField(obj, fid);
                    }

                    /**
                     * Allows to get the Java object field value
                     * @param env the JNI environment
                     * @param obj the object itself
                     * @param fid_name the field name
                     * @return the field value
                     */
                    bool get_bool_field(JNIEnv * env, jobject obj,
                            const char * fid_name) const {
                        jclass cls = env->GetObjectClass(obj);
                        jfieldID fid = env->GetFieldID(cls, fid_name, "Z");
                        return env->GetBooleanField(obj, fid);
                    }
                };
            }
        }
    }
}

#endif /* CONFIG_WRAPPER_HH */

