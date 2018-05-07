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
                typedef pair<double, double> min_max_pair;
                typedef vector<min_max_pair> min_max_cont;

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
                class scaled_ctrl_wrapper {
                private:
                    //Stores the data source reference
                    const data_source & m_data;
                    //Stores the controller class wrapper
                    ctrl_wrapper & m_wrapper;
                    //Stores the state index
                    int m_state_idx;
                    //Stores the computed scales per dof
                    double * m_scale;
                    //Stores the computed shifts per dof
                    double * m_shift;

                    /**
                     * Evaluates the controller on all the points and computes
                     * its minimum and maximum values.
                     * @param wrap the controller wrapper
                     * @param min_max the minimum/maximum pair per dof container
                     */
                    void eval_and_min_max_values(
                            ctrl_wrapper & wrap,
                            min_max_cont & min_max_data) {
                        //Default initialization for the min and max values
                        min_max_data.resize(m_wrapper.get_is_dim(),{DBL_MAX, 0.0});

                        //Instantiate the plain data iterator
                        sample_data_iter data_iter(m_data);

                        //Iterate over the available states
                        const double * state = NULL;
                        double input[m_wrapper.get_is_dim()];
                        while (data_iter.next_state(state)) {
                            //Compute the candidate's input for the state
                            wrap.compute_input(state, input);

                            //Update the min/max values if the input is good
                            for (int idx = 0; idx < m_wrapper.get_is_dim(); ++idx) {
                                min_max_data[idx].first =
                                        min(min_max_data[idx].first, input[idx]);
                                min_max_data[idx].second =
                                        max(min_max_data[idx].second, input[idx]);
                            }

                            //Skip on the inputs
                            data_iter.skip_inputs_data();
                        }
                    }

                public:

                    /**
                     * The basic constructor
                     * @param data the data source object
                     * @param wrap the regular wrapper
                     * @param scale the scale array pointer, NOT NULL
                     * @param shift the shift array pointer, NOT NULL
                     */
                    scaled_ctrl_wrapper(const data_source & data,
                            ctrl_wrapper & wrap,
                            double * scale, double * shift)
                    : m_data(data), m_wrapper(wrap), m_state_idx(0),
                    m_scale(scale), m_shift(shift) {
                        //Compute the inputs for the given individual and
                        //the min/max values per input dimension thereof
                        min_max_cont min_max;
                        eval_and_min_max_values(m_wrapper, min_max);

                        //Compute the scale and shift factors per dimension
                        for (int idx = 0; idx < wrap.get_is_dim(); ++idx) {
                            //Get the controller's and individuals min/max for the dof
                            const min_max_pair & orig_mm = m_data.m_min_max[idx];
                            const min_max_pair & ind_mm = min_max[idx];

                            //Compute the scaling factor
                            const double orig_span = (orig_mm.second - orig_mm.first);
                            const double int_span = (ind_mm.second - ind_mm.first);

                            //Set the scaling and shifting values
                            m_scale[idx] = (int_span == 0.0) ? 1.0 : orig_span / int_span;
                            m_shift[idx] = orig_mm.first - ind_mm.first * m_scale[idx];
                        }
                    }

                    /**
                     * The basic destructor
                     */
                    virtual ~scaled_ctrl_wrapper() {
                    }

                    /**
                     * When run in the caching mode must be reset to start
                     * evaluating states from the beginning.
                     */
                    inline void reset() {
                        m_state_idx = 0;
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
                     * @param input the input container to store the values into
                     */
                    inline void compute_input(const double *state, double *input) {
                        m_wrapper.compute_input(state, input);
                        for (int idx = 0; idx < m_wrapper.get_is_dim(); ++idx) {
                            input[idx] = m_scale[idx] * input[idx] + m_shift[idx];
                        }
                    }
                };
            }
        }
    }
}

#endif /* SCALED_CTRL_WRAPPER_HH */

