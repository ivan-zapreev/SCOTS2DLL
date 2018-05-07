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
#include <unordered_map>

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
                //Is used to store the set of distinct input values
                typedef unordered_set<uint64_t> inputs_set;
                //The vector of the input's abstract dof indexes
                typedef vector<uint64_t> abs_data;

                /**
                 * This class allows to perform the variance-based bisection of the hypercube.
                 * 1. Stores unique input id for counting inputs
                 * 2. Stores pointer to the abstract vector input
                 * 3. Computes input variance per component
                 * 4. Takes the norm of the variance vector
                 * 5. Uses this scalar value to choose bisection
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

                    //Stores the bisection right interval left border values
                    //per state-space dof or 0 if bisection is not possible
                    abs_data m_split_rlb;
                    //Stores the split sets per dof to count 
                    //inputs to choose the best bisection
                    vector<split_info> m_splits;

                    //Stores the concrete state classification in which
                    //bisection interval it passes, per dof
                    vector<split_side> m_state_split;
                    //This vector is used to store the distinct inputs 
                    //per dof per left/right bisection interval
                    vector<inputs_set> m_state_lr_ips;

                    //Is used to store the mapping of the input id into the input vector abstract states
                    unordered_map<uint64_t, abs_data> m_inputs;
                    //Stores the number of input space dimension
                    const int m_is_dim;

                    //The dof will be chosen after
                    int m_ch_dof_idx;
                    //The left/right sample size ratio
                    double m_lr_ratio;
                    //Stores the number of samples in the left
                    //hypercube for the chosen bisection
                    int64_t m_left_size;
                    //Stores the number of samples in the right
                    //hypercube for the chosen bisection
                    int64_t m_right_size;

                    /**
                     * Allows to split the interval in two sub-intervals
                     * @param lb the left interval border
                     * @param rb the right interval border
                     * @param lrb the right border of the new left interval
                     * @param rlb the left border of the new right interval
                     * @return true if the interval could be split, otherwise false
                     */
                    inline static bool bisect(
                            const uint64_t lb, const uint64_t rb,
                            uint64_t & lrb, uint64_t & rlb) {
                        bool result = true;

                        LOG("Bisecting interval [" << lb << ", " << rb << "]");

                        //If the number of intervals is less than 3 then we can not split
                        const uint64_t num_its = (rb - lb);
                        if (num_its >= 3) {
                            //Compute the mid point, is always rounded down
                            const uint64_t mid_point = (rb + lb) / 2;

                            LOG("Mid point of [" << lb << ", " << rb << "] is: " << mid_point);

                            //Compute the two new interval bounds
                            lrb = mid_point, rlb = mid_point;
                            //Check on the number of intervals
                            if (num_its % 2) {
                                //The number of intervals is odd, so the mid point
                                //is rounded down from the center, thus the left
                                //border of the right interval is to be on the
                                //right of the center
                                rlb += 1;
                            } else {
                                //The number of intervals is even, so there are
                                //two possibilities to make the bisection, so
                                //make a uniform random choice between them.
                                if (bool_gen()) {
                                    lrb -= 1;
                                } else {
                                    rlb += 1;
                                }
                            }
                            LOG("Splitting [" << lb << ", " << rb
                                    << "] into: [" << lb << ", " << lrb
                                    << "] and [" << rlb << ", " << rb << "]");
                        } else {
                            result = false;
                        }
                        return result;
                    }

                    /**
                     * Allows to increment the input count for the given input
                     * if it is seen for the first time for the given state.
                     * @param inputs the set of inputs for the given state
                     * @param cnt_map the map storing counts per input
                     * @param input the input value
                     */
                    void increment_count_if_new(
                            inputs_set & inputs, ips_counts & cnt_map, const uint64_t input) {
                        //Insert the input into the state inputs tracking set
                        const auto reg = inputs.insert(input);
                        //If this input has not been seen yet (the insertion was done)
                        if (reg.second) {
                            //Update the global count for the given input
                            auto it = cnt_map.find(input);
                            if (it == cnt_map.end()) {
                                cnt_map.insert({input, 1});
                            } else {
                                ++(it->second);
                            }
                        }
                    }

                    /**
                     * Allows to compute the root of a variance based on the given sample set.
                     * Computes the vector of variance values for each input vector component.
                     * Then a simple Euclidean norm of the vector is taken as the input's variance.
                     * The cubic root of the variance is taken as described in Chapter 7
                     * "Recursive Stratified Sampling" of the book:
                     * "Numerical Recipes in C: The Art of Scientific Computing".
                     * @param counts the input counts to compute the root variance from
                     * @param root_var the root of the variance to be set
                     * @return true if the variance could be computed, otherwise false
                     */
                    inline bool compute_root_variance(const ips_counts & counts, double & root_var) {
                        //See if we can compute the mean first
                        bool is_can_compute = (counts.size() > 0);
                        LOG("Inputs count: " << counts.size()
                                << ", compute mean (?): " << is_can_compute);

                        if (is_can_compute) {
                            //Compute the mean values vector
                            double sample_size = 0.0;
                            vector<double> mean_sum(m_is_dim, 0.0);
                            for (auto elem : counts) {
                                //Get the actual inputs vector from the input id
                                const abs_data & data = m_inputs[elem.first];
                                //Compute the mean sum as a vector, per coordinate
                                for (int idx = 0; idx < m_is_dim; ++idx) {
                                    //We take into account the number of times this input was seen 
                                    mean_sum[idx] += data[idx] * elem.second;
                                }
                                sample_size += elem.second;
                            }

                            //Check if the sample variance can be computed
                            is_can_compute = (sample_size > 1.0);
                            LOG("Sample count: " << sample_size
                                    << ", compute variance (?): " << is_can_compute);

                            if (is_can_compute) {
                                //Compute the partial sample variance
                                vector<double> var_sum(m_is_dim, 0.0);
                                for (auto elem : counts) {
                                    //Get the actual inputs vector from the input id
                                    const abs_data & data = m_inputs[elem.first];
                                    //Compute the variance sum as a vector, per coordinate
                                    for (int idx = 0; idx < m_is_dim; ++idx) {
                                        const double mean = mean_sum[idx] / sample_size;
                                        //We take into account the number of times this input was seen 
                                        var_sum[idx] += pow(data[idx] - mean, 2) * elem.second;
                                    }
                                }

                                //Compute the Euclidean norm of the variance vector
                                double norm_var = 0.0;
                                for (int idx = 0; idx < m_is_dim; ++idx) {
                                    const double var = var_sum[idx] / (sample_size - 1.0);
                                    norm_var += pow(var, 2);
                                }
                                norm_var = sqrt(norm_var);

                                //Compute the root of the variance
                                root_var = pow(norm_var, 1.0 / (1.0 + RSS_ALPHA));
                            }
                        }

                        return is_can_compute;
                    }

                public:

                    /**
                     * The basic constructor allowing to initialize the bisector based on the hypercube
                     * @param cube the hypercube to base the initialization on
                     * @param is_dim the number of input space dimensions
                     */
                    var_bisector(const random_hypercube & cube, const int is_dim)
                    : m_cube(cube), m_lr_ratio(NO_LR_RATIO), m_inputs(), m_is_dim(is_dim),
                    m_ch_dof_idx(UNDEFINED_DOF_IDX), m_left_size(0), m_right_size(0) {
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
                            if (bisect(ll[dof_idx], ur[dof_idx], lrb, rlb)) {
                                m_split_rlb[dof_idx] = rlb;
                            } else {
                                m_split_rlb[dof_idx] = NO_SPLIT_POINT;
                            }

                            LOG("The bisection left border of the right "
                                    << "cube is: " << m_split_rlb[dof_idx]);

                            //Now initialize the data
                            m_splits[dof_idx] = {
                                {0,
                                    {}},
                                {0,
                                    {}}
                            };
                            m_state_split[dof_idx] = split_side::UNDEF_SIDE;
                            m_state_lr_ips[dof_idx] = {};
                        }
                    }

                    /**
                     * Allows to start a new sample
                     * @param state_abs the abstract state
                     */
                    inline void start_sample(const abs_data & state_abs) {
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
                     * Allows to add a new sample state input.
                     * Stores unique input id for counting inputs
                     * Stores the corresponding abstract vector input values
                     * @param input the unique input id
                     * @param data the input's dof abstract coordinates
                     */
                    inline void add_input(uint64_t input, const abs_data & data) {
                        //Insert the input if it is not present yet
                        if (m_inputs.find(input) == m_inputs.end()) {
                            m_inputs[input] = data;
                        }

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
                        LOG("Finalizing variance bisector sampling.");
                        //Iterate over the dofs and choose the best one for bisection
                        for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx) {
                            LOG("Dof: " << dof_idx << " split left/right: "
                                    << m_split_rlb[dof_idx]);

                            //Check if there is a bisection possible in this dof
                            if (m_split_rlb[dof_idx] != NO_SPLIT_POINT) {
                                //Compute the bisection ratio
                                const split_info & lr_sets = m_splits[dof_idx];
                                const cube_info & left_cube = lr_sets.first;
                                const cube_info & right_cube = lr_sets.second;
                                double left_var = 0.0, right_var = 0.0;
                                //If both variance values can be computed
                                if (compute_root_variance(left_cube.second, left_var)
                                        && compute_root_variance(right_cube.second, right_var)) {
                                    const double lr_ratio = left_var / (left_var + right_var);
                                    //Check if the new ratio is larger than the old one
                                    if (abs(lr_ratio - NO_LR_RATIO) > abs(m_lr_ratio - NO_LR_RATIO)) {
                                        //Set the given dof as the currently chosen one
                                        m_lr_ratio = lr_ratio;
                                        m_ch_dof_idx = dof_idx;
                                        m_left_size = left_cube.first;
                                        m_right_size = right_cube.first;
                                    }
                                }
                            }
                        }
                    }

                    /**
                     * Allows to get the right and left borders of the left and
                     * right intervals for bisection for the chosen dof.
                     * @param lr the right border of the left interval to be set
                     * @param rl the left border of the right interval to be set
                     */
                    inline void get_lr_rl_dounds(uint64_t & lr, uint64_t & rl) {
                        rl = m_split_rlb[m_ch_dof_idx];
                        lr = rl - 1;
                    }

                    /**
                     * Must be called to indicate that all the samples have been
                     * added and one needs to decide in which dimension the
                     * bisection is to be done.
                     * @return the dimension in which bisection is to be done or
                     * negative value if such a dimension could not be found
                     */
                    inline int get_bisect_dof() {
                        return m_ch_dof_idx;
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

