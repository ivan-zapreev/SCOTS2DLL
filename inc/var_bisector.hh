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

#include <unordered_map>
#include <vector>
#include <cmath>
#include <unordered_set>

#include "jni_throw.hh"
#include "info_logger.hh"
#include "random_hypercube.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                typedef enum {
                    UNDEF_SIDE = 0,
                    LEFT_SIDE = UNDEF_SIDE + 1,
                    RIGHT_SIDE = LEFT_SIDE + 1,
                } split_side;
                //This data type matches the input abstract id in some input dof
                //into the number of times this id was observed in the sampled
                //input vectors, the input counts once per sampled state
                typedef unordered_map<uint64_t, uint64_t> ips_counts;
                //Stores number of samples in the hypercube mapped to the corresponding input counts
                typedef pair<int64_t, ips_counts> cube_info;
                //Stores input info for the left and right cubes for the split
                typedef pair<cube_info, cube_info> split_info;
                typedef unordered_set<uint64_t> inputs_set;

                /**
                 * This class allows to perform the variance-based bisection of the hypercube
                 */
                class var_bisector {
                public:
                    //Defines the value for the undefined dof index
                    static constexpr int UNDEFINED_DOF_IDX = -1;
                    //The value for the left bound of the right interval for when no split is to be done in the dof
                    static constexpr int64_t NO_SPLIT_POINT = 0;
                    //Defines the RSS method lambda parameter for sample ratio computation
                    static constexpr double RSS_ALPHA = 2.0;
                    //Defines the value for the left-right sample ratio for which no bisection is to be done
                    static constexpr double NO_LR_RATIO = 0.5;

                private:
                    //Stores the reference to the corresponding random hypercube
                    const random_hypercube & m_cube;
                    //Stores the input space dof index for 
                    //which the variance bisection will be done
                    const int m_is_dof_idx;

                    //Stores the bisection right interval left border values
                    //per state-space dof or 0 if bisection is not possible
                    vector<uint64_t> m_split_rlb;
                    //Stores the split sets per dof to count 
                    //inputs to choose the best bisection
                    vector<split_info> m_splits;

                    //Stores the concrete state classification in which
                    //bisection interval it passes, per dof
                    vector<split_side> m_state_split;
                    //This vector is used to store the distinct inputs 
                    //per dof per left/right bisection interval
                    vector<inputs_set> m_state_lr_ips;

                    //The dof will be chosen after
                    int m_chosen_dof;
                    double m_lr_ratio;
                    int64_t m_left_size;
                    int64_t m_right_size;

                    /**
                     * Allows to increment the input count for the given set
                     * @param set the set to work with
                     * @param input the input whoes count is to be incremented
                     */
                    void increment_count_if_new(inputs_set & inputs, ips_counts & set, const uint64_t input) {
                        //Insert the input into the state inputs tracking set
                        const auto reg = inputs.insert(input);
                        //If this input has not been seen yet (the insertion was done)
                        if (reg.second) {
                            //Update the global inputs count
                            auto it = set.find(input);
                            if (it == set.end()) {
                                set.insert({input, 1});
                            } else {
                                ++(it->second);
                            }
                        }
                    }

                    /**
                     * Allows to compute the root of variance based on the given sample set
                     * @param counts the input counts to compute the root variance from
                     * @param root_var the root of the variance to be set
                     * @return true if the variance could be computed, otherwise false
                     */
                    inline bool get_root_variance(const ips_counts & counts, double & root_var) {
                        //Check if the sample variance can be computed
                        const bool is_can_compute = (counts.size() >= 2);

                        //If it can be computed then do the computations
                        if (is_can_compute) {
                            //Compute the sample mean value
                            double count = 0.0;
                            double mean = 0.0;
                            for (auto elem : counts) {
                                mean += elem.first * elem.second;
                                count += elem.second;
                            }
                            mean = mean / count;

                            //Compute the sample variance
                            root_var = 0.0;
                            for (auto elem : counts) {
                                const double delta = (elem.first * elem.second - mean);
                                root_var += delta * delta;
                            }
                            root_var = root_var / (count - 1.0);

                            //Compute the root of the variance
                            root_var = pow(root_var, 1.0 / (1.0 + RSS_ALPHA));
                        }

                        return is_can_compute;
                    }

                public:

                    /**
                     * The basic constructor allowing to initialize the bisector based on the hypercube
                     * @param cube the hypercube to base the initialization on
                     * @param is_dof_idx the input dof index for variance computations
                     */
                    var_bisector(const random_hypercube & cube, const int is_dof_idx)
                    : m_cube(cube), m_is_dof_idx(is_dof_idx), m_lr_ratio(NO_LR_RATIO),
                    m_chosen_dof(UNDEFINED_DOF_IDX), m_left_size(0), m_right_size(0) {
                        //Get the number of dofs
                        const int num_dofs = m_cube.get_num_dim();
                        LOG("Constructing a new variance bisector for " << num_dofs << " dofs");

                        //Reserve the container elements
                        m_split_rlb.resize(num_dofs);
                        m_splits.resize(num_dofs);
                        m_state_split.resize(num_dofs);
                        m_state_lr_ips.resize(num_dofs);

                        //Initialize the split left interval right borders
                        uint64_t lrb = 0, rlb = 0;
                        const uint64_t * ll = cube.get_ll();
                        const uint64_t * ur = cube.get_ur();
                        for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx) {
                            LOG("Considering bisection dof: " << dof_idx);
                            
                            //Initialization with 0 of right interval left border
                            //means there is no possible split in this dimension
                            if (random_hypercube::split_interval(
                                    ll[dof_idx], ur[dof_idx], lrb, rlb)) {
                                m_split_rlb[dof_idx] = rlb;
                            } else {
                                m_split_rlb[dof_idx] = NO_SPLIT_POINT;
                            }
                            
                            LOG("The bisection left border of the right "
                                    << "cube is: " << m_split_rlb[dof_idx]);

                            //Now initialize the data
                            m_split_rlb[dof_idx] = NO_SPLIT_POINT;
                            m_splits[dof_idx] = {{0,{}},{0,{}}};
                            m_state_split[dof_idx] = split_side::UNDEF_SIDE;
                            m_state_lr_ips[dof_idx] = {};
                        }
                    }

                    /**
                     * Allows to start a new sample
                     * @param state_abs the abstract state
                     */
                    inline void start_sample(const vector<uint64_t> & state_abs) {
                        //Get the number of dofs
                        const int num_dofs = m_cube.get_num_dim();
                        //Iterate over all dofs and decide into which side of
                        //the split the input space dof value will be put
                        for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx) {
                            if (m_split_rlb[dof_idx] != NO_SPLIT_POINT) {
                                if (state_abs[dof_idx] >= m_split_rlb[dof_idx]) {
                                    m_state_split[dof_idx] = split_side::RIGHT_SIDE;
                                } else {
                                    m_state_split[dof_idx] = split_side::LEFT_SIDE;
                                }
                            } else {
                                m_state_split[dof_idx] = split_side::UNDEF_SIDE;
                            }
                        }
                    }

                    /**
                     * Allows to add a new sample state input
                     * @param input_abs the abstract input
                     */
                    inline void add_input(const vector<uint64_t> & input_abs) {
                        //Get the input to be classified
                        const uint64_t input = input_abs[m_is_dof_idx];
                        //Get the number of dofs
                        const int num_dofs = m_cube.get_num_dim();
                        //Iterate over the dofs and register it in a proper sets
                        for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx) {
                            inputs_set & inputs = m_state_lr_ips[dof_idx];
                            split_info & lr_sets = m_splits[dof_idx];
                            //Increment the element counts if it is a new inputs for the sample
                            switch (m_state_split[dof_idx]) {
                                case split_side::LEFT_SIDE:
                                    increment_count_if_new(
                                            inputs, lr_sets.first.second, input);
                                    break;
                                case split_side::RIGHT_SIDE:
                                    increment_count_if_new(
                                            inputs, lr_sets.second.second, input);
                                    break;
                                default:
                                    //Nothing to be done the split in this dimension is not possible
                                    break;
                            }
                        }
                    }

                    /**
                     * Allows to indicate the end of the sample
                     */
                    inline void stop_sample() {
                        //Get the number of dofs
                        const int num_dofs = m_cube.get_num_dim();
                        for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx) {
                            //For the state there must have been always an input,
                            //so we need add one to the left or right sample count
                            split_info & lr_sets = m_splits[dof_idx];
                            switch (m_split_rlb[dof_idx]) {
                                case split_side::LEFT_SIDE:
                                    //Increment the left sample count for the dof 
                                    ++lr_sets.first.first;
                                    break;
                                case split_side::RIGHT_SIDE:
                                    //Increment the right sample count for the dof 
                                    ++lr_sets.second.first;
                                    break;
                                default:
                                    //Nothing to be done as there was not split in this dof
                                    break;
                            }
                            //Clear the inputs set as it will be re-used
                            m_state_lr_ips[dof_idx].clear();
                        }
                    }

                    /**
                     * Will be called once the sampling for the
                     * given bisector is finished.
                     */
                    inline void done_sampling() {
                        //Get the number of dofs
                        const int num_dofs = m_cube.get_num_dim();
                        //Iterate over the dofs and choose the best one for bisection
                        for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx) {
                            //Check if there is a bisection possible in this dof
                            if (m_split_rlb[dof_idx] != NO_SPLIT_POINT) {
                                //Compute the bisection ratio
                                const split_info & lr_sets = m_splits[dof_idx];
                                const cube_info & left_cube = lr_sets.first;
                                const cube_info & right_cube = lr_sets.second;
                                double left_var = 0.0, right_var = 0.0;
                                //If both variance values can be computed
                                if (get_root_variance(left_cube.second, left_var)
                                        && get_root_variance(right_cube.second, right_var)) {
                                    const double lr_ratio = left_var / (left_var + right_var);
                                    //Check if the new ratio is larger than the old one
                                    if (abs(lr_ratio - NO_LR_RATIO) > abs(m_lr_ratio - NO_LR_RATIO)) {
                                        //Set the given dof as the currently chosen one
                                        m_lr_ratio = lr_ratio;
                                        m_chosen_dof = dof_idx;
                                        m_left_size = left_cube.first;
                                        m_right_size = right_cube.first;
                                    }
                                }
                            }
                        }
                    }

                    /**
                     * Must be called to indicate that all the samples have been
                     * added and one needs to decide in which dimension the
                     * bisection is to be done.
                     * @return the dimension in which bisection is to be done or
                     * negative value if such a dimension could not be found
                     */
                    inline int choose_dof() {
                        return m_chosen_dof;
                    }

                    /**
                     * Get the variance-based ratio for deciding which part of 
                     * the samples to be samples from the left sub-cube
                     * @return a value from [0.0..1.0] which gives the proportion
                     * of samples to use for the left sub-cube
                     */
                    inline double get_lr_ratio() {
                        return m_lr_ratio;
                    }

                    /**
                     * Allows to get the number of samples in the left hypercube
                     * @return the number of samples in the left hypercube
                     */
                    inline int64_t get_left_size() {
                        return m_left_size;
                    }

                    /**
                     * Allows to get the number of samples in the right hypercube
                     * @return the number of samples in the right hypercube
                     */
                    inline int64_t get_right_size() {
                        return m_right_size;
                    }
                };
            }
        }
    }
}


#endif /* VAR_BISECTOR_HH */

