/* 
 * File:   bdd_data_iter.hh
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
 * Created on April 27, 2018, 3:35 PM
 */

#ifndef BDD_DATA_ITER_HH
#define BDD_DATA_ITER_HH

#include <vector>

#include "jni_throw.hh"
#include "info_logger.hh"
#include "config_wrapper.hh"
#include "data_source.hh"
#include "data_iter.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * This class represents the bdd data iterator used in Monte Carlo simulations
                 */
                class bdd_data_iter : public data_iter {
                private:
                    //Stores the configuration reference
                    const config_wrapper & m_config;
                    //The dummy inputs array
                    raw_data m_dummy;
                    //Will store the state inputs
                    raw_data const * m_p_inputs;
                    //Store the current inputs iterator
                    raw_data::const_iterator m_curr_input;
                    //The vector to store a single input in the original space
                    raw_data m_input;
                    //The array to store a single input in the discrete space
                    double * m_ainput;

                public:

                    /**
                     * The basic constructor
                     * 
                     * @param data the data source to be used
                     */
                    bdd_data_iter(const data_source & data)
                    : data_iter(data), m_config(data.get_config()), m_dummy(),
                    m_p_inputs(&m_dummy), m_curr_input(m_dummy.end()),
                    m_input(m_config.m_num_is_dim),
                    m_ainput(new double[m_config.m_num_is_dim]) {
                    }

                    /**
                     * The basic destructor
                     */
                    virtual ~bdd_data_iter() {
                        if (m_ainput) {
                            delete[] m_ainput;
                        }
                    }

                    /**
                     * Allows to pre-load inputs for the given state. 
                     * @param state the state to pre-load inputs for.
                     * @return true if the inputs were pre-loaded, otherwise false
                     */
                    bool extract_inputs(const uint64_t id) {
                        //Get the list of inputs
                        m_p_inputs = &m_data.restriction(id);

                        //Get the begin and end iterators
                        m_curr_input = m_p_inputs->begin();

                        //Return false if there is not inputs, otherwise true
                        return (m_p_inputs->size() != 0);
                    }

                    /**
                     * Allows to get the next abstract input (grid indexes) from
                     * the container in a form of a vector of double values.
                     * 
                     * @param input the reference to the container pointer to be set
                     * @return true if the next input exists, otherwise false
                     */
                    inline bool next_input(const double * & input) {
                        const bool has_input = (m_curr_input != m_p_inputs->end());
                        if (has_input) {
                            //Get a new input
                            m_input.assign(m_curr_input, m_curr_input + m_config.m_num_is_dim);

                            //Convert state to abstract state and add to data
                            m_data.m_p_is_mgr->xtois(m_input, m_ainput);

                            //Set the input data pointer
                            input = m_ainput;

                            //Move forward in the list of inputs
                            m_curr_input += m_config.m_num_is_dim;
                        } else {
                            input = NULL;
                        }
                        return has_input;
                    }

                    /**
                     * Once the state is read we can skip to the next state by
                     * skipping the inputs corresponding to the last read states.
                     */
                    inline void skip_inputs_data() {
                        //Just set the current input to be the end input
                        m_curr_input = m_p_inputs->end();
                    }
                };
            }
        }
    }
}

#endif /* INPUTS_ITER_HH */

