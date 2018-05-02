/* 
 * File:   var_bisector.hh
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
 * Created on May 1, 2018, 4:35 PM
 */

#ifndef VAR_BISECTOR_HH
#define VAR_BISECTOR_HH

#include <vector>

#include "jni_throw.hh"
#include "info_logger.hh"
#include "random_hypercube.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * This class allows to perform the variance-based bisection of the hypercube
                 */
                class var_bisector {
                public:

                    /**
                     * The basic constructor allowing to initialize the bisector based on the hypercube
                     * @param cube the hypercube to base the initialization on
                     * @param is_dof_idx the input dof index for variance computations
                     */
                    var_bisector(const random_hypercube & cube, const int is_dof_idx) {
                        //ToDo: Implement
                    }

                    /**
                     * Allows to start a new sample
                     * @param state_abs the abstract state
                     */
                    void start_sample(const vector<uint64_t> & state_abs) {
                        //ToDo: Implement
                    }

                    /**
                     * Allows to add a new sample state input
                     * @param input_abs the abstract input
                     */
                    void add_input(const vector<uint64_t> & input_abs) {
                        //ToDo: Implement
                    }

                    /**
                     * Allows to finish adding the sample
                     */
                    void finish_sample() {
                        //ToDo: Implement
                    }

                    /**
                     * Must be called to indicate that all the samples have been
                     * added and one needs to decide in which dimension the
                     * bisection is to be done.
                     * @return the dimension in which bisection is to be done or
                     * negative value if such a dimension could not be found
                     */
                    int choose_dof() {
                        //ToDo: Implement
                        return -1;
                    }

                    /**
                     * Get the variance-based ratio for deciding which part of 
                     * the samples to be samples from the left sub-cube
                     * @return a value from [0.0..1.0] which gives the proportion
                     * of samples to use for the left sub-cube
                     */
                    double get_lvr_ratio() {
                        //ToDo: Implement
                        return 0.5;
                    }

                    /**
                     * Allows to get the number of samples in the left hypercube
                     * @return the number of samples in the left hypercube
                     */
                    int64_t get_left_size() {
                        //ToDo: Implement
                        return 0;
                    }

                    /**
                     * Allows to get the number of samples in the right hypercube
                     * @return the number of samples in the right hypercube
                     */
                    int64_t get_right_size() {
                        //ToDo: Implement
                        return 0;
                    }
                };
            }
        }
    }
}


#endif /* VAR_BISECTOR_HH */

