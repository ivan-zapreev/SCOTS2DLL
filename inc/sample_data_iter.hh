/* 
 * File:   plain_data_iter.hh
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
 * Created on April 26, 2018, 9:09 AM
 */

#ifndef SAMPLE_DATA_ITER_HH
#define SAMPLE_DATA_ITER_HH

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
                 * This class allows to iterate through the controller's sample data.
                 */
                class sample_data_iter : public data_iter {
                private:
                    //Stores the current cursor
                    double const * m_cursor;
                    //Stores the current state inputs count
                    abs_type m_ips_cnt;
                    //Stores the current state input index
                    abs_type m_ips_idx;

                    /**
                     * When we enter the inputs part, re-initialize input count and index
                     */
                    inline void prepare_for_inputs() {
                        //Get the inputs count
                        m_ips_cnt = static_cast<abs_type> (*m_cursor);
                        //Set the inputs index to the first one
                        m_ips_idx = 0;
                        //Shift to the inputs part
                        m_cursor += 1;
                    }

                public:

                    /**
                     * The basic constructor
                     * 
                     * @param data the data source to be used
                     */
                    sample_data_iter(const data_source & data)
                    : data_iter(data) {
                        reset();
                    }

                    /**
                     * Allows to retrieve the sample size
                     * @return the sample size
                     */
                    inline abs_type get_sample_size() const {
                        return m_data.get_sample_size();
                    }

                    /**
                     * Allows to reset the iterator to start from the beginning
                     */
                    inline void reset() {
                        m_cursor = m_data.get_sample_start();
                        m_ips_cnt = 0;
                        m_ips_idx = 0;
                    }

                    /**
                     * Allows to check if we have iterated through all the data yet, or not
                     * 
                     * @return true if there is still data to iterate
                     */
                    inline bool has_data() {
                        return (m_cursor < m_data.get_sample_end());
                    }

                    /**
                     * Allows to get the next state from the container.
                     * 
                     * @param state the reference to the container pointer to be set
                     * @return true if the next state exists, otherwise false
                     */
                    inline bool next_state(const double * & state) {
                        bool has_next = has_data();
                        if (has_next) {
                            //Get the current cursor position
                            state = m_cursor;
                            //Shift the cursor to where the inputs count is
                            m_cursor += m_config.m_num_ss_dim;
                            //Prepare for inputs iteration
                            prepare_for_inputs();
                        }
                        return has_next;
                    }

                    /**
                     * Allows to get the next input from the container.
                     * 
                     * @param input the reference to the container pointer to be set
                     * @return true if the next input exists, otherwise false
                     */
                    inline bool next_input(const double * & input) {
                        bool has_next = (m_ips_idx < m_ips_cnt);
                        //Iterate to the next input if it is present
                        if (has_next) {
                            //Set the input
                            input = m_cursor;
                            //Increment the index
                            m_ips_idx++;
                            //Shift to the following input
                            m_cursor += m_config.m_num_is_dim;
                        }
                        return has_next;
                    }

                    /**
                     * Allows to skip the state.
                     * 
                     * @return true if the state existed, otherwise false
                     */
                    inline void skip_state_data() {
                        //Shift the cursor to where the inputs count is
                        m_cursor += m_config.m_num_ss_dim;
                        //Prepare for inputs iteration
                        prepare_for_inputs();
                    }

                    /**
                     * Once the state is read we can skip to the next state by
                     * skipping the inputs corresponding to the last read states.
                     * This method also works if some of the inputs have already been read.
                     */
                    inline void skip_inputs_data() {
                        //Shift through the inputs part
                        m_cursor += m_config.m_num_is_dim * (m_ips_cnt - m_ips_idx);
                    }

                };
            }
        }
    }
}
#endif /* SAMPLE_DATA_ITER_HH */

