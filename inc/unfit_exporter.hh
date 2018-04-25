/* 
 * File:   unfit_exporter.hh
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
 * Created on April 24, 2018, 8:26 AM
 */

#ifndef UNFIT_EXPORTER_HH
#define UNFIT_EXPORTER_HH


#include "info_logger.hh"
#include "config_wrapper.hh"
#include "data_source.hh"
#include "ctrl_wrapper.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * This class is used for exporting the unfit points
                 */
                class unfit_exporter : public error_computer {
                private:
                    //Stores the unfit states BDD of the controller
                    BDD m_unfit_bdd;

                public:

                    unfit_exporter(info_logger & log, data_source & data)
                    : error_computer(log, data), m_unfit_bdd() {
                    }

                    /**
                     * The basic destructor
                     */
                    virtual ~unfit_exporter() {
                        //De-couple the BDD from the cudd manager
                        m_unfit_bdd = BDD();
                    }

                    /**
                     * Allows to start the unfit points export by re-setting
                     * the internal data storage for the unfit points
                     * @param env the environment
                     */
                    void start_unfit_export(JNIEnv * env) {
                        LOG(m_log, "Starting unfit points export");
                        //Re-initialize the unfit points BDD
                        m_unfit_bdd = m_data.m_p_cudd->bddZero();
                        //Make sure all the data points are pre-loaded,
                        //this is needed to speed up computations
                        m_data.ensure_data_points();
                    }

                    /**
                     * Allows to add the unfit points of this individual into
                     * the unfit points to be exported into the BDD
                     * @param env the JNI environment
                     * @param wrap the controller wrapper
                     * @return the real (actual) fitness of the individual
                     */
                    double compute_unfit_points(JNIEnv * env, ctrl_wrapper & wrap) {
                        LOG(m_log, "Starting unfit points export of " << wrap.get_name()
                                << " for dof " << wrap.get_dof_idx());
                        double result = 0.0;

                        if (m_config.m_num_ss_dim > 0) {
                            const int dof_idx = wrap.get_dof_idx();
                            if ((dof_idx >= 0) && (dof_idx < (m_config.m_num_is_dim))) {
                                result = add_unfit_points_to_bdd(env, wrap);
                            } else {
                                (void) throwException(m_log, env, IllegalArgumentException,
                                        "Improper input space dimension index!");
                            }
                        } else {
                            (void) throwException(m_log, env, IllegalStateException,
                                    "The state-space size is not set!");
                        }

                        return result;
                    }

                    /**
                     * Allows to store the collected so-far unfit points into the BDD.
                     * @param env the JNI environment
                     * @param file_name the BDD file name
                     */
                    void finish_unfit_export(JNIEnv * env, const char * file_name) {
                        LOG(m_log, "Dumping unfit points into " << file_name);
                        store_controller(*m_data.m_p_cudd,
                                m_data.m_p_ss_mgr->get_states_set(),
                                m_unfit_bdd, string(file_name));
                        //De-couple the BDD from the CUDD manager
                        m_unfit_bdd = BDD();
                    }

                private:

                    /**
                     * Allows to store the unfit points into a BDD.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     */
                    double add_unfit_points_to_bdd(JNIEnv * env, ctrl_wrapper & wrap) {
                        //Get the dof index relative to the largest state dof index
                        const int is_dof_idx = wrap.get_dof_idx();
                        //The variables to store the number of fit and unfit points
                        long total = 0.0, unfit = 0.0;

                        //Define the cursor
                        double const * cursor = &m_data.m_all_data[0];
                        abs_type state_abs[m_config.m_num_ss_dim];
                        while (cursor < m_data.m_cursor_max) {
                            //Compute the state error
                            double delta_err = DBL_MAX;
                            const double ctr_value = wrap.compute_input(cursor, m_config.m_num_ss_dim);
                            double const * state_dbl = cursor;
                            //If there are nun of inf values or the error delta is
                            //larger than the maximum then this point is not fit.
                            if (!compute_state_error(ctr_value, is_dof_idx, cursor, delta_err)
                                    || (delta_err >= MAX_ATTRACTOR_SIZE)) {
                                //Copy the values into the state array
                                for (int idx = 0; idx < m_config.m_num_ss_dim; ++idx) {
                                    state_abs[idx] = *state_dbl;
                                    ++state_dbl;
                                }

                                //Convert the abstract state to BDD and add to the total
                                m_unfit_bdd = m_unfit_bdd | m_data.m_p_ss_mgr->x_to_bdd(state_abs);

                                //Increment the number of unfit points
                                ++unfit;
                            }
                            //Increment the number of fit points
                            ++total;
                        }

                        LOG(m_log, "The number of unfit points for dof "
                                << wrap.get_dof_idx() << " is " << unfit
                                << "/" << total);

                        //Return the actual fitness value
                        return ((double) (total - unfit)) / ((double) total);
                    }

                };

            }
        }
    }
}

#endif /* UNFIT_EXPORTER_HH */

