/* 
 * File:   data_manager.hh
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
 * Created on April 23, 2018, 10:42 AM
 */

#ifndef DATA_SOURCE_HH
#define DATA_SOURCE_HH

#define SCOTS_BDD 1
#include "scots.hh"

#include "inputs_mgr.hh"
#include "states_mgr.hh"
#include "input_output.hh"
#include "jni_throw.hh"

#include "info_logger.hh"
#include "config_wrapper.hh"

using namespace std;
using namespace scots;
using namespace tud::ctrl::scots::optimal;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                typedef vector<double> raw_data;
                typedef vector<abs_type> abs_data;

                /**
                 * This class is used to load and store the controller's data
                 */
                class data_source {
                public:
                    //Stores the CUDD manager
                    Cudd * m_p_cudd;
                    //Stores the BDD of the controller
                    BDD * m_p_ctrl_bdd;
                    //Stores the Symbolic set of the controller
                    SymbolicSet * m_p_ctr;

                    //Stores the inputs manager
                    inputs_mgr * m_p_is_mgr;
                    //Stores the states manager
                    states_mgr * m_p_ss_mgr;

                    //Stores the number of the grid data points
                    abs_type m_num_states;
                    //The continuous array of states
                    vector<double> m_all_data;
                    //Stores the maximum cursor value
                    double * m_cursor_max;

                    //Stores the min/max values for the input dofs
                    vector<pair<double, double>> m_min_max;

                    /**
                     * The basic constructor.
                     * @param log the logger object to be used
                     */
                    data_source(info_logger & log)
                    : m_log(log), m_p_cudd(NULL), m_p_ctrl_bdd(NULL),
                    m_p_ctr(NULL), m_p_is_mgr(NULL), m_p_ss_mgr(NULL),
                    m_num_states(0), m_all_data(),
                    m_cursor_max(NULL), m_min_max() {
                    }

                    virtual ~data_source() {
                        if (m_p_ctr) {
                            delete m_p_ctr;
                        }
                        if (m_p_ctrl_bdd) {
                            delete m_p_ctrl_bdd;
                        }
                        if (m_p_ss_mgr) {
                            delete m_p_ss_mgr;
                        }
                        if (m_p_is_mgr) {
                            delete m_p_is_mgr;
                        }
                        if (m_p_cudd) {
                            delete m_p_cudd;
                        }
                    }

                    /**
                     * Loads the controller with the given name.
                     * Must be called once per instance
                     * @param env the JNI environment
                     * @param file_name the file name
                     * @return the controller dimensions
                     */
                    int load(JNIEnv * const env, const char * file_name) {
                        //Read controller from file
                        LOG(m_log, "Start initializing data structures");

                        //Make a new BDD
                        if (m_p_ctrl_bdd) delete m_p_ctrl_bdd;
                        m_p_ctrl_bdd = new BDD();

                        //Make a new symbolic set
                        if (m_p_ctr) delete m_p_ctr;
                        m_p_ctr = new SymbolicSet();

                        //Create a new CUDD manager
                        if (m_p_cudd) delete m_p_cudd;
                        m_p_cudd = new Cudd();

                        //Disable automatic variable ordering
                        m_p_cudd->AutodynDisable();

                        //Read controller from file
                        LOG(m_log, "Start loading the controller: " << file_name);

                        //Load the controller
                        bool is_fail = false;
                        if (!read_from_file(*m_p_cudd, *m_p_ctr, *m_p_ctrl_bdd, file_name)) {
                            is_fail = true;
                        }
                        LOG(m_log, "Loading: " << file_name << " is finished");

                        //Throw if could not load otherwise extract the bdd data
                        if (is_fail) {
                            (void) throwException(m_log, env, FileNotFoundException,
                                    "Failed to read the controller file!");
                        }

                        return m_p_ctr->get_dim();
                    }

                    /**
                     * Allows to get the number of points of the state-space grid.
                     * @param env the JNI environment
                     * @param ss_dim the number of state-space dimensions
                     * @return the number of points of the state-space grid
                     */
                    int get_state_space_size(JNIEnv * env, jint ss_dim) {
                        LOG(m_log, "Retrieving the number of points on the state-space grid");
                        int size = 1;
                        if (m_p_ctr != NULL) {
                            const vector<abs_type> nu_gp = m_p_ctr->get_no_gp_per_dim();
                            for (int idx = 0; idx < ss_dim; ++idx) {
                                size *= nu_gp[idx];
                            }
                        } else {
                            (void) throwException(m_log, env, IllegalStateException,
                                    "The controller is not loaded yet!");
                        }
                        return size;
                    }

                    /**
                     * Allows to configure the fitness computer based on the configuration object
                     * @param env the JNI environment
                     * @param cfg the configuration object
                     */
                    void configure(JNIEnv * env, jobject cfg) {
                        //Create a new configuration wrapper
                        config_wrapper config(env, cfg, m_p_ctr->get_dim());

                        //Check that the configuration is correct
                        config.finalize(m_log, env);

                        //Store the reference to the configuration object
                        m_config = config;

                        //Clean the old data if any
                        re_initialize_data();

                        //Convert the points into the internal data structures, only if needed!
                        if (m_config.m_is_scale) {
                            if (m_config.m_is_mc) {
                                convert_points_to_data<true, false>();
                            } else {
                                convert_points_to_data<true, true>();
                            }
                        } else {
                            if (m_config.m_is_mc) {
                                convert_points_to_data<false, false>();
                            } else {
                                convert_points_to_data<false, true>();
                            }
                        }
                    }

                    /**
                     * Allows to retrieve the data source configuration.
                     * 
                     * @return the data source configuration
                     */
                    const config_wrapper & get_config() const {
                        return m_config;
                    }

                    /**
                     * Allows to ensure that the data points are loaded
                     */
                    void ensure_data_points() {
                        if (m_all_data.size() == 0) {
                            convert_points_to_data<false, true>();
                        }
                    }

                private:
                    //Stores the reference to the logging object
                    info_logger & m_log;

                    //Stores the configuration object
                    config_wrapper m_config;

                    void re_initialize_data() {
                        //Initialize the min-max pairs
                        m_min_max.clear();

                        //Clear the old data points
                        m_all_data.clear();

                        //Delete any previously present managers
                        if (m_p_ss_mgr) delete m_p_ss_mgr;
                        if (m_p_is_mgr) delete m_p_is_mgr;

                        //Initialize the input and state manager
                        m_p_is_mgr = new inputs_mgr(*m_p_ctr, m_config.m_num_ss_dim);
                        m_p_ss_mgr = new states_mgr(*m_p_ctr, m_config.m_num_ss_dim,
                                *m_p_ctrl_bdd, *m_p_cudd,
                                m_p_is_mgr->get_inputs_set());
                    }

                    /**
                     * Allows to convert the points array into the internal data structures
                     * @param num_pg the number of points in the array
                     * @param arr_gp the points array
                     * @return the number of states
                     */
                    template<bool IS_SCALE_FIT, bool IS_EXPLICIT_DATA>
                    void convert_points_to_data() {
                        LOG(m_log, "Start extracting grid points");

                        //Initialize the min/max values
                        if (IS_SCALE_FIT) {
                            for (abs_type idx = 0; idx < m_config.m_num_is_dim; ++idx) {
                                m_min_max.emplace_back(DBL_MAX, 0.0);
                            }
                        }

                        if (IS_SCALE_FIT || IS_EXPLICIT_DATA) {
                            //Get the number the states with inputs
                            const raw_data all_states = m_p_ss_mgr->get_points();

                            //Define state and abstract state containers
                            raw_data state(m_config.m_num_ss_dim), input(m_config.m_num_is_dim);
                            abs_data astate(m_config.m_num_ss_dim), ainput(m_config.m_num_is_dim);

                            //Iterate over the states, get the corresponding
                            //inputs and add them to the estimator set by ids
                            auto state_begin = all_states.begin();
                            while (state_begin != all_states.end()) {
                                //Get a new state vector
                                state.assign(state_begin, state_begin + m_config.m_num_ss_dim);

                                //Get the list of inputs
                                raw_data state_inputs = m_p_ctr->restriction(*m_p_cudd, *m_p_ctrl_bdd, state);

                                //If we need to store all data explicitly for numeric fitness computations
                                if (IS_EXPLICIT_DATA) {
                                    //Convert state to abstract state and add to data
                                    m_p_ss_mgr->xtois(state, astate);
                                    m_all_data.insert(m_all_data.end(), astate.begin(), astate.end());

                                    //Get the number of inputs and add to the data
                                    const int ips_cnt = state_inputs.size() / m_config.m_num_is_dim;
                                    m_all_data.push_back(ips_cnt);
                                }

                                //Iterate over the inputs
                                auto input_begin = state_inputs.begin();
                                while (input_begin != state_inputs.end()) {
                                    //Get a new input vector
                                    input.assign(input_begin, input_begin + m_config.m_num_is_dim);

                                    //Convert state to abstract state and add to data
                                    m_p_is_mgr->xtois(input, ainput);
                                    //If we need to store all data explicitly for numeric fitness computations
                                    if (IS_EXPLICIT_DATA) {
                                        m_all_data.insert(m_all_data.end(), ainput.begin(), ainput.end());
                                    }

                                    //Compute the min and max values per input dimension
                                    if (IS_SCALE_FIT) {
                                        for (abs_type idx = 0; idx < m_config.m_num_is_dim; ++idx) {
                                            pair<double, double> & elem = m_min_max[idx];
                                            elem.first = min(elem.first, static_cast<double> (ainput[idx]));
                                            elem.second = max(elem.second, static_cast<double> (ainput[idx]));
                                        }
                                    }

                                    //Move forward in the list of inputs
                                    input_begin += m_config.m_num_is_dim;
                                }

                                //Move forward in the list of states
                                state_begin += m_config.m_num_ss_dim;
                            }
                        }

                        //Add the MAX_STATE_DEVIATION to the bounds
                        if (IS_SCALE_FIT) {
                            for (abs_type idx = 0; idx < m_config.m_num_is_dim; ++idx) {
                                m_min_max[idx].first -= MAX_STATE_DEVIATION;
                                m_min_max[idx].second += MAX_STATE_DEVIATION;
                            }
                        }

                        //Set the maximum cursor value
                        if (IS_EXPLICIT_DATA) {
                            m_cursor_max = &m_all_data[m_all_data.size() - 1] + 1;
                        } else {
                            m_cursor_max = NULL;
                        }

                        //Set the number of states value
                        m_num_states = m_p_ss_mgr->get_size();
                        LOG(m_log, "The state-space size is = " << m_num_states);
                    }
                };
            }
        }
    }
}

#endif /* DATA_SOURCE_HH */

