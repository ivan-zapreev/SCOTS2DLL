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
#include "sample_data_iter.hh"

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

                    unfit_exporter(data_source & data)
                    : error_computer(data), m_unfit_bdd() {
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
                        LOG("Starting unfit points export");
                        //Re-initialize the unfit points BDD
                        m_unfit_bdd = m_data.m_p_cudd->bddZero();
                    }

                    /**
                     * Allows to add the unfit points of this individual into
                     * the unfit points to be exported into the BDD
                     * @param env the JNI environment
                     * @param wrap the controller wrapper
                     * @return the real (actual) fitness of the individual
                     */
                    double compute_unfit_points(JNIEnv * env, ctrl_wrapper & wrap) {
                        LOG("Starting unfit points export of " << wrap.get_name());
                        double result = 0.0;

                        if (m_config.m_num_ss_dim > 0) {
                            result = add_unfit_points_to_bdd(env, wrap);
                        } else {
                            (void) throwException(env, IllegalStateException,
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
                        LOG("Dumping unfit points into " << file_name);
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
                        //Get controller's states
                        raw_data states;
                        const abs_type dom_size = m_data.m_p_ss_mgr->get_size();
                        LOG("Start getting " << dom_size << " domain states");
                        m_data.m_p_ss_mgr->get_points(states);

                        //Convert the current individual into a BDD I
                        LOG("Start converting the symbolic controller into BDD");

                        auto state_begin = states.begin();
                        BDD ctrl_bdd = m_data.m_p_cudd->bddZero();

                        raw_data state(m_config.m_num_ss_dim);
                        raw_data astate(m_config.m_num_ss_dim);
                        raw_data ainput(m_config.m_num_is_dim);

                        int64_t new_cnt = 0, old_cnt = 0, state_cnt = 0;

                        //Iterate over the states, get the corresponding
                        while (state_begin != states.end()) {
                            //Get a new state vector
                            state.assign(state_begin, state_begin + m_config.m_num_ss_dim);

                            //Convert state to an abstract one
                            m_data.m_p_ss_mgr->xtois(state, astate);

                            //Compute the state input on the abstract state
                            wrap.compute_input(astate.data(), ainput.data());

                            //Add the state/input pair to the controller BDD
                            ctrl_bdd |=
                                    (m_data.m_p_ss_mgr->i_to_bdd(astate)
                                    & m_data.m_p_is_mgr->i_to_bdd(ainput));

                            //Move forward in the list of states
                            state_begin += m_config.m_num_ss_dim;

                            //De the logging for convenience
                            new_cnt = state_cnt * 1000 / dom_size;
                            if (new_cnt != old_cnt) {
                                LOG("Converted " << ((double) new_cnt) / 10.0
                                        << "% of the symbolic controller");
                                old_cnt = new_cnt;
                            }
                            ++state_cnt;
                        }

                        LOG("Start computing the good part of the symbolic controller");
                        //Get the controller BDD and intersect it with the original
                        BDD good_ctrl_bdd = ctrl_bdd & *m_data.m_p_ctrl_bdd;

                        //Get the domains of both the original and the 
                        //good part of the generated controller
                        BDD U = m_data.m_p_is_mgr->get_inputs_set().get_cube(*m_data.m_p_cudd);
                        LOG("Start computing the original domain");
                        BDD dom_orig = m_data.m_p_ctrl_bdd->ExistAbstract(U);
                        LOG("Start computing the good part domain");
                        BDD dom_good = good_ctrl_bdd.ExistAbstract(U);

                        //The unfit states are those which are not in the good part
                        LOG("Start computing the unfit domain part");
                        m_unfit_bdd = dom_orig & !dom_good;

                        LOG("Start computing the number of unfit points in BDD");
                        const abs_type unfit_size
                                = m_data.m_p_ss_mgr->get_states_set().get_size(
                                *m_data.m_p_cudd, m_unfit_bdd);

                        //Return the actual fitness value
                        LOG("The number of unfit points is " << unfit_size << "/" << dom_size);
                        return ((double) (dom_size - unfit_size)) / ((double) dom_size);
                    }
                };
            }
        }
    }
}

#endif /* UNFIT_EXPORTER_HH */

