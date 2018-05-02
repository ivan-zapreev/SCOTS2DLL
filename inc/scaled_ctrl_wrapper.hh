/* 
 * File:   scaled_ctrl_wrapper.hh
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
 * Created on April 26, 2018, 1:45 PM
 */

#ifndef SCALED_CTRL_WRAPPER_HH
#define SCALED_CTRL_WRAPPER_HH

#include "info_logger.hh"
#include "config_wrapper.hh"
#include "data_source.hh"
#include "sample_data_iter.hh"
#include "ctrl_wrapper.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * This class represents a generic scaled CTRL wrapper.
                 * In the caching mode requires more memory but is faster.
                 * In the non-scaling mode is more memory efficient but slower.
                 * 
                 * NOTE: In the caching mode the provided state argument is not 
                 * used and just the next pre-computed input is given for the 
                 * same order of states as are present in the plain data of the
                 * data_source. Once all states are evaluated one needs to reset
                 * the class in order to start evaluation from the start again.
                 * 
                 * @param IS_CACHING if true then the input values are 
                 * pre-computed while scaling factors are computed and
                 * are stored explicitly.
                 */
                template<bool IS_CACHING>
                class scaled_ctrl_wrapper {
                private:
                    //Stores the data source reference
                    const data_source & m_data;
                    //Stores the controller class wrapper
                    ctrl_wrapper & m_wrapper;
                    //Stores the state index
                    int m_state_idx;
                    //Stores the array of inputs
                    double * m_inputs;
                    //Stores the computed scale
                    double m_scale;
                    //Stores the computed shift
                    double m_shift;

                    /**
                     * Evaluates the controller on all the points and computes
                     * its minimum and maximum values.
                     * @param data_iter the data iterator
                     * @param wrap the controller wrapper
                     * @param min_input the minimum value to compute
                     * @param max_input the maximum value to compute
                     * @return a new array of point value or NULL if the 
                     * controller evaluation on some of the points was
                     * infinite or nan.
                     */
                    double * eval_and_min_max_values(
                            sample_data_iter & data_iter, ctrl_wrapper & wrap,
                            double &min_input, double &max_input) {
                        //Allocate the inputs array and initialize the state index
                        abs_type state_idx = 0;
                        double * inputs = (IS_CACHING ? new double[data_iter.get_sample_size()] : NULL);

                        //Default initialization for the min and max values
                        max_input = 0.0;
                        min_input = DBL_MAX;

                        //Iterate over the available states
                        const double * state = NULL;
                        while (data_iter.next_state(state)) {
                            //Compute the candidate's input for the state
                            const double ctr_value = wrap.compute_input(state);
                            //Check if the value is sensible
                            //Store the value
                            if (IS_CACHING) {
                                inputs[state_idx] = ctr_value;
                            }

                            //Update the min/max values if the input is good
                            if (!isnan(ctr_value) && !isinf(ctr_value)) {
                                max_input = max(max_input, ctr_value);
                                min_input = min(min_input, ctr_value);
                            }

                            //Skip on the inputs
                            data_iter.skip_inputs_data();
                            //Move to the next state index
                            ++state_idx;
                        }
                        return inputs;
                    }

                public:

                    /**
                     * The basic constructor
                     * @param data the data source object
                     * @param wrapper the regular wrapper
                     * @param scale the scaling factor
                     * @param shift the shifting factor
                     */
                    scaled_ctrl_wrapper(const data_source & data, ctrl_wrapper & wrapper)
                    : m_data(data), m_wrapper(wrapper), m_state_idx(0) {
                        //Get the input dof index
                        const int is_dof_idx = m_wrapper.get_dof_idx();

                        //Instantiate the plain data iterator
                        sample_data_iter data_iter(m_data, is_dof_idx);
                        
                        //Compute the inputs for the given individual and the min/max values thereof
                        double min_input, max_input;
                        m_inputs = eval_and_min_max_values(
                                data_iter, m_wrapper, min_input, max_input);

                        //Get the original controller min/max values
                        const pair<double, double> & min_max = m_data.m_min_max[is_dof_idx];
                        //Compute the scaling factor
                        const double orig_span = (min_max.second - min_max.first);
                        const double func_span = (max_input - min_input);

                        //Set the scaling and shifting values                        
                        m_scale = (func_span == 0.0) ? 1.0 : orig_span / func_span;
                        m_shift = min_max.first - min_input * m_scale;
                    }

                    /**
                     * The basic destructor
                     */
                    virtual ~scaled_ctrl_wrapper() {
                        if (IS_CACHING && m_inputs) {
                            delete[] m_inputs;
                        }
                    }

                    /**
                     * When run in the caching mode must be reset to start
                     * evaluating states from the beginning.
                     */
                    inline void reset() {
                        m_state_idx = 0;
                    }

                    /**
                     * Allows to retrieve the computed scaling factor
                     * @return  the computed scaling factor
                     */
                    inline double get_scale() const {
                        return m_scale;
                    }

                    /**
                     * Allows to retrieve the computed shifting factor
                     * @return  the computed shifting factor
                     */
                    inline double get_shift() const {
                        return m_shift;
                    }

                    /**
                     * Allows to get the input dof index
                     * @return the input dof index
                     */
                    inline const int & get_dof_idx() const {
                        return m_wrapper.get_dof_idx();
                    }

                    /**
                     * Allows to get the wrapper's name
                     * @return the wrapper's name
                     */
                    inline const string & get_name() const {
                        return m_wrapper.get_name();
                    }

                    /**
                     * Computes the input value for the given state
                     * @param state the state to call the method for
                     * @return the computed input value for this state
                     */
                    inline double compute_input(const double *state) {
                        double ctrl_val = (IS_CACHING ? m_inputs[m_state_idx++] : m_wrapper.compute_input(state));
                        return ctrl_val * m_scale + m_shift;
                    }
                };
            }
        }
    }
}

#endif /* SCALED_CTRL_WRAPPER_HH */

