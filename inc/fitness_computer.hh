/*
 * File:   ftn_computer.hh
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

#define PI 3.14159265

#include "info_logger.hh"
#include "config_wrapper.hh"
#include "data_source.hh"
#include "ctrl_wrapper.hh"
#include "error_computer.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                typedef priority_queue<double, vector<double>, std::greater<double>> scores_queue;

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
                    template<bool is_complex, fitness_type ftn_type>
                    void update_scores(JNIEnv * env, const double delta_err,
                            double & exact_score, double & part_score,
                            scores_queue & part_scores) {
                        exact_score += ((delta_err < MAX_ATTRACTOR_SIZE) ? 1.0 : 0.0);
                        if (is_complex) {
                            if (delta_err <= m_config.m_attr_size) {
                                part_score += 1.0;
                            } else {
                                switch (ftn_type) {
                                    case EXACT_FITNESS:
                                        break;
                                    case ATANG_FITNESS:
                                        part_scores.push(1.0 - 2.0 * atan((delta_err - m_config.m_attr_size) * m_config.m_ftn_scale) / PI);
                                        break;
                                    case INVER_FITNESS:
                                        part_scores.push(1.0 / (1.0 + (delta_err - m_config.m_attr_size) * m_config.m_ftn_scale));
                                        break;
                                    default:
                                        (void) throwException(m_log, env,
                                                IllegalArgumentException,
                                                "Unsupported fitness type!");
                                }
                            }
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
                    template<bool is_complex, fitness_type ftn_type>
                    void compute_fitness_scale(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        //Allocate the inputs array
                        double * inputs = new double[m_data.m_num_states];
                        //Set the cursor
                        double const * cursor = &m_data.m_all_data[0];
                        //Iterate over controller states and compute inputs
                        abs_type state_idx = 0;
                        double max_input = 0.0;
                        double min_input = DBL_MAX;
                        while (cursor < m_data.m_cursor_max) {
                            //Compute the candidate's input for the state
                            const double ctr_value = wrap.compute_input(cursor, m_config.m_num_ss_dim);
                            //Check if the value is sensible
                            if (!isnan(ctr_value) && !isinf(ctr_value)) {
                                //Store the value
                                inputs[state_idx] = ctr_value;
                                //Compute the min/max values
                                max_input = max(max_input, ctr_value);
                                min_input = min(min_input, ctr_value);
                                //Shift the cursor to where the inputs count is
                                cursor += m_config.m_num_ss_dim;
                                //Get the inputs count
                                const abs_type ips_cnt = static_cast<abs_type> (*cursor);
                                //Shift through the inputs part
                                cursor += 1 + m_config.m_num_is_dim * ips_cnt;
                                //Move to the next state index
                                ++state_idx;
                            } else {
                                //The value is off so it makes no sense to proceed
                                LOG(m_log, "WARNING: Dof " << wrap.get_dof_idx()
                                        << " individual: " << wrap.get_name()
                                        << " input value is: " << ctr_value);
                                delete[] inputs;
                                return;
                            }
                        }

                        //Get the dof index relative to the largest state dof index
                        const int is_dof_idx = wrap.get_dof_idx();
                        //Initialize the priority queue for safe computations
                        scores_queue part_scores;
                        //Get the original controller min/max values
                        const pair<double, double> & min_max = m_data.m_min_max[is_dof_idx];
                        //Compute the scaling factor
                        const double orig_span = (min_max.second - min_max.first);
                        const double func_span = (max_input - min_input);

                        //Set the scaling and shifting values                        
                        scale = (func_span == 0.0) ? 1.0 : orig_span / func_span;
                        shift = min_max.first - min_input * scale;

                        //Re-initialize the cursor
                        cursor = &m_data.m_all_data[0];
                        //Initialize the scores
                        double part_score_sum = 0.0;
                        double exact_score = 0.0;
                        double part_score = 0.0;
                        //Iterate over actual controller and compute errors
                        state_idx = 0;
                        while (cursor < m_data.m_cursor_max) {
                            //Compute the state error
                            double delta_err = DBL_MAX;
                            const double adj_input = inputs[state_idx] * scale + shift;
                            if (compute_state_error(adj_input, is_dof_idx, cursor, delta_err)) {
                                //Compute the scores
                                update_scores<is_complex, ftn_type>(env, delta_err,
                                        exact_score, part_score, part_scores);
                            } else {
                                //We still got an error that is off, so we give up 
                                delete[] inputs;
                                return;
                            }
                            ++state_idx;
                        }

                        //Compute the resulting fitness values
                        while (!part_scores.empty()) {
                            part_score_sum += part_scores.top();
                            part_scores.pop();
                        }

                        exact_ftn = exact_score / m_data.m_num_states;
                        if (is_complex) {
                            const double ext_ftn = (part_score + part_score_sum) / m_data.m_num_states;
                            req_ftn = sqrt(exact_ftn * exact_ftn + ext_ftn * ext_ftn);
                        } else {
                            req_ftn = exact_ftn;
                        }
                        //Clear memory
                        delete[] inputs;
                    }

                    /**
                     * Allows to compute the error of the given controller.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     */
                    template<bool is_complex, fitness_type ftn_type>
                    void compute_fitness_plain(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn) {
                        //Get the dof index relative to the largest state dof index
                        const int is_dof_idx = wrap.get_dof_idx();
                        //Re-set the values to zero first
                        double exact_score = 0.0;
                        double part_score = 0.0;
                        //Create an inverted priority queue
                        scores_queue part_scores;

                        //Define the cursor
                        double const * cursor = &m_data.m_all_data[0];
                        double part_score_sum = 0.0;
                        while (cursor < m_data.m_cursor_max) {
                            //Compute the state error
                            double delta_err = DBL_MAX;
                            const double ctr_value = wrap.compute_input(cursor, m_config.m_num_ss_dim);
                            if (compute_state_error(ctr_value, is_dof_idx, cursor, delta_err)) {
                                //Compute the scores
                                update_scores<is_complex, ftn_type>(env, delta_err,
                                        exact_score, part_score, part_scores);
                            } else {
                                //We still got an error that is off, so we give up 
                                return;
                            }
                        }

                        //Compute the resulting fitness values
                        while (!part_scores.empty()) {
                            part_score_sum += part_scores.top();
                            part_scores.pop();
                        }

                        exact_ftn = exact_score / m_data.m_num_states;
                        if (is_complex) {
                            const double ext_ftn = (part_score + part_score_sum) / m_data.m_num_states;
                            req_ftn = sqrt(exact_ftn * exact_ftn + ext_ftn * ext_ftn);
                        } else {
                            req_ftn = exact_ftn;
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
                    void compute_fitness(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        switch (m_config.m_ftn_type) {
                            case EXACT_FITNESS:
                                if (m_config.m_is_scale) {
                                    compute_fitness_scale< IS_COMPLEX, EXACT_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn, scale, shift);
                                } else {
                                    compute_fitness_plain< IS_COMPLEX, EXACT_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn);
                                }
                                break;
                            case ATANG_FITNESS:
                                if (m_config.m_is_scale) {
                                    compute_fitness_scale< IS_COMPLEX, ATANG_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn, scale, shift);
                                } else {
                                    compute_fitness_plain< IS_COMPLEX, ATANG_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn);
                                }
                                break;
                            case INVER_FITNESS:
                                if (m_config.m_is_scale) {
                                    compute_fitness_scale< IS_COMPLEX, INVER_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn, scale, shift);
                                } else {
                                    compute_fitness_plain< IS_COMPLEX, INVER_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn);
                                }
                                break;
                            default:
                                (void) throwException(m_log, env,
                                        IllegalArgumentException,
                                        "Unsupported fitness type!");
                        }
                    }

                public:

                    /**
                     * The basic constructor
                     * 
                     * @param log the reference to the logging object
                     * @param data the configured data source
                     */
                    fitness_computer(info_logger & log, data_source & data) :
                    error_computer(log, data) {
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
                                    compute_fitness<true>(env, wrap, exact_ftn, req_ftn, scale, shift);
                                } else {
                                    compute_fitness<false>(env, wrap, exact_ftn, req_ftn, scale, shift);
                                }
                            } else {
                                (void) throwException(m_log, env, IllegalArgumentException,
                                        "Improper input space dimension index!");
                            }
                        } else {
                            (void) throwException(m_log, env, IllegalStateException,
                                    "The state-space size is not set!");
                        }
                    }
                };
            }
        }
    }
}

#endif /* FITNESS_COMPUTER_HH */

