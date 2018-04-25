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

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * The base class for fitness error computing entities
                 */
                class error_computer {
                protected:
                    //Stores the reference to the logging object
                    info_logger & m_log;

                    //Stores the data source reference
                    data_source & m_data;

                    //Stores the configuration object reference
                    const config_wrapper & m_config;

                    /**
                     * Allows to compute the state error given the input from the candidate
                     * @param ctr_value the control value for the state
                     * @param is_dof_idx the dof index in the dof dimensions
                     * @param cursor the cursor for the original controller data
                     * @param delta_err the error to be computed
                     * @return if true then the error could be computed, otherwise false
                     */
                    bool compute_state_error(const double ctr_value, const int is_dof_idx,
                            double const * & cursor, double & delta_err) {
                        bool result = true;

                        //Do not proceed if there is a nan value
                        if (!isnan(ctr_value) && !isinf(ctr_value)) {
                            //Shift the cursor to where the inputs count is
                            cursor += m_config.m_num_ss_dim;
                            //Get the inputs count
                            const abs_type ips_cnt = static_cast<abs_type> (*cursor);
                            //Shift to the inputs part
                            cursor += 1;
                            //Now get the minimum error
                            for (abs_type idx = 0; idx < ips_cnt; ++idx) {
                                //Get the dof input value 
                                const double act_value = cursor[is_dof_idx];
                                //Compute the delta error
                                delta_err = min(delta_err, abs(ctr_value - act_value));
                                //Shift to the following input
                                cursor += m_config.m_num_is_dim;
                            }
                        } else {
                            LOG(m_log, "ERROR: Individual control value is: " << ctr_value);
                            result = false;
                        }
                        return result;
                    }

                public:

                    /**
                     * The basic constructor
                     * 
                     * @param log the reference to the logging object
                     * @param data the configured data source
                     */
                    error_computer(info_logger & log, data_source & data)
                    : m_log(log), m_data(data), m_config(m_data.get_config()) {
                    }

                };

            }
        }
    }
}

#endif /* ERROR_COMPUTER_HH */
