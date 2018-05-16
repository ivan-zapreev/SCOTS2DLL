/*
 * Copyright (C) 2018 Dr. Ivan S. Zapreev <ivan.zapreev@gmail.com>
 *
 *  Visit my Linked-in profile:
 *     https://nl.linkedin.com/in/zapreevis
 *  Visit my GitHub:
 *     https://github.com/ivan-zapreev
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * File:   error_computer.hh
 * Author: Dr. Ivan S. Zapreev <ivan.zapreev@gmail.com>
 *
 * Created on April 24, 2018, 8:42 AM
 */

#ifndef ERROR_COMPUTER_HH
#define ERROR_COMPUTER_HH

#include "info_logger.hh"
#include "config_wrapper.hh"
#include "data_source.hh"
#include "error_computer.hh"
#include "sample_data_iter.hh"

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * The base class for fitness error computing entities
                 */
                class error_computer {
                protected:
                    //Stores the data source reference
                    data_source & m_data;

                    //Stores the configuration object reference
                    const config_wrapper & m_config;

                    //Stores the number of input space dimensions
                    const int m_is_dim;

                    /**
                     * Allows to compute the state error given the input from the candidate
                     * @param input_cnt the control input vector values for the state as computed by the individual
                     * @param ipts_iter the inputs iterator for the original controller data
                     * @param delta_err the error to be computed
                     * @return if true then the error could be computed, otherwise false
                     */
                    bool min_state_error(const double * ind_input,
                            sample_data_iter & ipts_iter, double & delta_err) {
                        //Initialize the error
                        delta_err = DBL_MAX;

                        //Check that the vector has only proper inputs
                        for (int idx = 0; idx < m_is_dim; ++idx) {
                            //If the individual's input is not valid then
                            if (isnan(ind_input[idx]) || isinf(ind_input[idx])) {
                                //Skip the inputs data
                                ipts_iter.skip_inputs_data();
                                //And return false
                                return false;
                            }
                        }

                        //If the input is good then compute its error
                        //Iterate over the inputs and compute the minimum 
                        //(over all controller's state inputs) maximum 
                        //(over all input vector dimensions) delta error
                        const double * act_input = NULL;
                        while (ipts_iter.next_input(act_input)) {
                            //Compute the maximum dof error
                            double max_dof_err = 0.0;
                            for (int idx = 0; idx < m_is_dim; ++idx) {
                                const double dof_error = abs(ind_input[idx] - act_input[idx]);
                                max_dof_err = max(max_dof_err, dof_error);
                            }
                            //Compute the delta error
                            delta_err = min(delta_err, max_dof_err);
                        }
                        return true;
                    }

                public:

                    /**
                     * The basic constructor
                     * 
                     * @param data the configured data source
                     */
                    error_computer(data_source & data)
                    : m_data(data), m_config(m_data.get_config()),
                    m_is_dim(m_config.m_num_is_dim) {
                    }

                };

            }
        }
    }
}

#endif /* ERROR_COMPUTER_HH */

