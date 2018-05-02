/*
 * File:   fitness_computer.hh
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
 * Created on January 14, 2018, 8:32 PM
 */

#ifndef FITNESS_COMPUTER_HH
#define FITNESS_COMPUTER_HH

#include <cmath>
#include <string>
#include <queue>
#include <functional>
#include <cstring>
#include <limits>
#include <vector>
#include <utility>
#include <ctime>
#include <unordered_set>

#define PI 3.14159265

#include "info_logger.hh"
#include "config_wrapper.hh"
#include "data_source.hh"
#include "ctrl_wrapper.hh"
#include "scaled_ctrl_wrapper.hh"
#include "error_computer.hh"
#include "random_hypercube.hh"
#include "sample_data_iter.hh"
#include "bdd_data_iter.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                typedef priority_queue<double, vector<double>, std::greater<double>> scores_queue;
                typedef unordered_set<uint64_t> sample_set;

                /**
                 * The fitness compute class
                 */
                class fitness_computer : public error_computer {
                protected:

                    /**
                     * Allows to update the scores given the delta error
                     * @param env the jni environment
                     * @param delta_err the delta error
                     * @param exact_score the exact score to update
                     * @param part_score the integer part of the partial score
                     * @param part_scores the decimal scores queue
                     */
                    template<bool IS_COMPLEX, fitness_type FTN_TYPE>
                    void update_scores(JNIEnv * env, const double delta_err,
                            double & exact_score, double & part_score,
                            scores_queue & part_scores) {
                        exact_score += ((delta_err < MAX_ATTRACTOR_SIZE) ? 1.0 : 0.0);
                        if (IS_COMPLEX) {
                            if (delta_err <= m_config.m_attr_size) {
                                part_score += 1.0;
                            } else {
                                switch (FTN_TYPE) {
                                    case EXACT_FITNESS:
                                        break;
                                    case ATANG_FITNESS:
                                        part_scores.push(1.0 - 2.0 * atan((delta_err - m_config.m_attr_size) * m_config.m_ftn_scale) / PI);
                                        break;
                                    case INVER_FITNESS:
                                        part_scores.push(1.0 / (1.0 + (delta_err - m_config.m_attr_size) * m_config.m_ftn_scale));
                                        break;
                                    default:
                                        (void) throwException(env,
                                                IllegalArgumentException,
                                                "Unsupported fitness type!");
                                }
                            }
                        }
                    }

                    /**
                     * Allows to compute finalize the fitness computations.
                     * @param exact_score the exact score value
                     * @param part_score the partial score value
                     * @param part_scores the partial scores
                     * @param num_points the number of points
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     */
                    template<bool IS_COMPLEX>
                    void finalize_fitness_computations(
                            const double & exact_score, const double & part_score,
                            scores_queue & part_scores, const abs_type num_points,
                            double & exact_ftn, double & req_ftn) {
                        //Check if the number of points is not zero
                        if (num_points > 0) {
                            //Compute the resulting fitness values
                            double part_score_sum = 0.0;
                            while (!part_scores.empty()) {
                                part_score_sum += part_scores.top();
                                part_scores.pop();
                            }

                            exact_ftn = exact_score / num_points;
                            if (IS_COMPLEX) {
                                const double ext_ftn = (part_score + part_score_sum) / num_points;
                                req_ftn = sqrt(exact_ftn * exact_ftn + ext_ftn * ext_ftn);
                            } else {
                                req_ftn = exact_ftn;
                            }
                        } else {
                            LOG("WARNING: The number of computed fitness points is 0!");
                        }
                    }

                    /**
                     * Allows to compute the error of the given controller 
                     * using numerical methods and without scaling.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     */
                    template<bool IS_COMPLEX, fitness_type FTN_TYPE,
                    typename WRAP_TYPE, typename IPTS_ITER_TYPE>
                    void compute_min_state_error(
                            JNIEnv * env, IPTS_ITER_TYPE & data_iter,
                            WRAP_TYPE & wrap, const int is_dof_idx,
                            const double * state, double & delta_err,
                            double & exact_score, double & part_score,
                            scores_queue & part_scores) {
                        const double ctr_value = wrap.compute_input(state);
                        if (min_state_error<IPTS_ITER_TYPE>(ctr_value, is_dof_idx, data_iter, delta_err)) {
                            //Compute the scores
                            update_scores<IS_COMPLEX, FTN_TYPE>(env, delta_err,
                                    exact_score, part_score, part_scores);
                        }
                    }

                    /**
                     * Allows to compute the error of the given controller 
                     * using numerical methods and without scaling.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     */
                    template<bool IS_COMPLEX, fitness_type FTN_TYPE, typename WRAP_TYPE>
                    void compute_fitness_pl(JNIEnv * env, WRAP_TYPE & wrap,
                            double & exact_ftn, double & req_ftn) {
                        //Get the dof index relative to the largest state dof index
                        const int is_dof_idx = wrap.get_dof_idx();

                        //Re-set the values to zero first
                        double exact_score = 0.0;
                        double part_score = 0.0;

                        //Create an inverted priority queue
                        scores_queue part_scores;

                        //Instantiate the plain data iterator
                        sample_data_iter data_iter(m_data, is_dof_idx);

                        //Iterate over the stages and compute errors
                        const double * state = NULL;
                        double delta_err = DBL_MAX;
                        while (data_iter.next_state(state)) {
                            //Compute the state error
                            compute_min_state_error<IS_COMPLEX, FTN_TYPE>(
                                    env, data_iter, wrap, is_dof_idx, state,
                                    delta_err, exact_score, part_score, part_scores);
                        }

                        //Finalize the fitness computations
                        finalize_fitness_computations<IS_COMPLEX>(
                                exact_score, part_score, part_scores,
                                m_data.get_sample_size(is_dof_idx),
                                exact_ftn, req_ftn);
                    }

                    /**
                     * Allows to compute the error of the given controller
                     * using numerical methods with scaling.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     * @param scale scale for the case scaling is requested
                     * @param shift shift for the case scaling is requested
                     */
                    template<bool IS_COMPLEX, fitness_type FTN_TYPE>
                    void compute_fitness_sc(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        //Initialize the scaled wrapper
                        scaled_ctrl_wrapper<false> scaled_wrap(m_data, wrap);

                        //Compute the fitness
                        compute_fitness_pl<IS_COMPLEX, FTN_TYPE>(
                                env, scaled_wrap, exact_ftn, req_ftn);

                        //Get the scaling values
                        scale = scaled_wrap.get_scale();
                        shift = scaled_wrap.get_shift();
                    }

                    /**
                     * Allows to compute the error of the given controller.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     * @param scale scale for the case scaling is requested
                     * @param shift shift for the case scaling is requested
                     */
                    template<bool IS_COMPLEX, bool IS_SCALE, fitness_type FTN_TYPE>
                    void compute_cmpl_sc_ft(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        if (IS_SCALE) {
                            compute_fitness_sc<IS_COMPLEX, FTN_TYPE>(
                                    env, wrap, exact_ftn, req_ftn, scale, shift);
                        } else {
                            compute_fitness_pl<IS_COMPLEX, FTN_TYPE>(
                                    env, wrap, exact_ftn, req_ftn);
                        }
                    }

                    /**
                     * Allows to compute the error of the given controller.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     * @param scale scale for the case scaling is requested
                     * @param shift shift for the case scaling is requested
                     */
                    template<bool IS_COMPLEX, bool IS_SCALE>
                    void compute_cmpl_sc(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        switch (m_config.m_ftn_type) {
                            case EXACT_FITNESS:
                                compute_cmpl_sc_ft< IS_COMPLEX, IS_SCALE, EXACT_FITNESS >
                                        (env, wrap, exact_ftn, req_ftn, scale, shift);
                                break;
                            case ATANG_FITNESS:
                                compute_cmpl_sc_ft< IS_COMPLEX, IS_SCALE, ATANG_FITNESS >
                                        (env, wrap, exact_ftn, req_ftn, scale, shift);
                                break;
                            case INVER_FITNESS:
                                compute_cmpl_sc_ft< IS_COMPLEX, IS_SCALE, INVER_FITNESS >
                                        (env, wrap, exact_ftn, req_ftn, scale, shift);
                                break;
                            default:
                                (void) throwException(env,
                                        IllegalArgumentException,
                                        "Unsupported fitness type!");
                        }
                    }

                    /**
                     * Allows to compute the error of the given controller.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     * @param scale scale for the case scaling is requested
                     * @param shift shift for the case scaling is requested
                     */
                    template<bool IS_COMPLEX>
                    void compute_cmpl(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        if (m_config.m_is_scale) {
                            compute_cmpl_sc< IS_COMPLEX, true >
                                    (env, wrap, exact_ftn, req_ftn, scale, shift);
                        } else {
                            compute_cmpl_sc< IS_COMPLEX, false >
                                    (env, wrap, exact_ftn, req_ftn, scale, shift);
                        }
                    }

                public:

                    /**
                     * The basic constructor
                     * 
                     * @param data the configured data source
                     */
                    fitness_computer(data_source & data) :
                    error_computer(data) {
                    }

                    /**
                     * The basic destructor
                     */
                    virtual ~fitness_computer() {
                    }

                    /**
                     * Allows to compute fitness of the given controller wrapper
                     * @param env the JNI environment
                     * @param wrap the controller wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     * @param scale scale for the case scaling is requested
                     * @param shift shift for the case scaling is requested
                     */
                    void compute(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        if (m_config.m_num_ss_dim > 0) {
                            const int dof_idx = wrap.get_dof_idx();
                            if ((dof_idx >= 0) && (dof_idx < m_config.m_num_is_dim)) {
                                if (m_config.m_is_complex) {
                                    compute_cmpl<true>(env, wrap, exact_ftn, req_ftn, scale, shift);
                                } else {
                                    compute_cmpl<false>(env, wrap, exact_ftn, req_ftn, scale, shift);
                                }
                            } else {
                                (void) throwException(env, IllegalArgumentException,
                                        "Improper input space dimension index!");
                            }
                        } else {
                            (void) throwException(env, IllegalStateException,
                                    "The state-space size is not set!");
                        }
                    }
                };
            }
        }
    }
}

#endif /* FITNESS_COMPUTER_HH */

