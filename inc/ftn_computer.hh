/*
 * File:   ctrl_wrapper.hh
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
 * Created on January 14, 2018, 8:32 PM
 */

#ifndef FTN_COMPUTER_HH
#define FTN_COMPUTER_HH

#include <cmath>
#include <string>
#include <iostream>
#include <queue>
#include <functional>
#include <cstring>
#include <limits>
#include <vector>
#include <utility>
#include <ctime>

#define PI 3.14159265

#define SCOTS_BDD 1
#include "scots.hh"

#include "inputs_mgr.hh"
#include "states_mgr.hh"
#include "input_output.hh"
#include "jni_throw.hh"
#include "ctrl_wrapper.hh"

using namespace std;
using namespace scots;
using namespace tud::ctrl::scots::optimal;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                typedef priority_queue<double, vector<double>, std::greater<double>> scores_queue;
                typedef vector<double> raw_data;
                typedef vector<abs_type> abs_data;

                /**
                 * The fitness types enumeration
                 */
                enum fitness_type {
                    UNDEF_FITNESS = 0,
                    EXACT_FITNESS = 1,
                    ATANG_FITNESS = 2,
                    INVER_FITNESS = 3
                };

                /**
                 * The fitness compute class
                 */
                class ftn_computer {
                private:
                    //Stores the maximum attractor size
                    static constexpr double MAX_ATTRACTOR_SIZE = 0.5;
                    static constexpr double EPSILON = 1e-2;
                    //Stores the maximum stare deviation
                    static constexpr double MAX_STATE_DEVIATION = (MAX_ATTRACTOR_SIZE - EPSILON);

                    //Stores the log file
                    ofstream m_log_file;
                    //Stores the CUDD manager
                    Cudd * m_p_cudd;
                    //Stores the BDD of the controller
                    BDD * m_p_ctrl_bdd;
                    //Stores the inputs manager
                    inputs_mgr * m_p_is_mgr;
                    //Stores the states manager
                    states_mgr * m_p_ss_mgr;
                    //Stores the Symbolic set of the controller
                    SymbolicSet * m_p_ctr;
                    //Stores the unfit states BDD of the controller
                    BDD m_unfit_bdd;
                    //Stores the state-space size
                    int m_ss_dim;
                    //Stores the input-space size
                    int m_is_dim;
                    //Stores the number of the grid data points
                    abs_type m_num_states;
                    //The continuous array of states
                    vector<double> m_all_data;
                    //Stores the maximum cursor value
                    double * m_cursor_max;
                    //Stores the attractor size
                    double m_attr_size;
                    //Stores the min/max values for the input dofs
                    vector<pair<double, double>> m_min_max;
                    //Stores the scaling flag
                    bool m_is_scale;
                    //Stores the fitness type
                    fitness_type m_ftn_type;
                    //Stores the flag for complex fitness
                    bool m_is_complex;
                    //Stores the fitness function scaling factor for the complex fitness types INVERSE and ATANGENT
                    double m_ftn_scale;

                    //The logging operator
                    friend ostream& operator<<(ftn_computer& strm, const char * str);

                protected:

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
                            cursor += m_ss_dim;
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
                                cursor += m_is_dim;
                            }
                        } else {
                            *this << "ERROR: Individual control value is: " << ctr_value;
                            result = false;
                        }
                        return result;
                    }

                    /**
                     * Allows to update the scores given the delta error
                     * @param env the jni environment
                     * @param delta_err the delta error
                     * @param exact_score the exact score to update
                     * @param part_score the integer part of the partial score
                     * @param part_scores the decimal scores queue
                     */
                    template<bool is_complex, fitness_type ftn_type>
                    void update_scores(JNIEnv * env, const double delta_err,
                            double & exact_score, double & part_score,
                            scores_queue & part_scores) {
                        exact_score += ((delta_err < MAX_ATTRACTOR_SIZE) ? 1.0 : 0.0);
                        if (is_complex) {
                            if (delta_err <= m_attr_size) {
                                part_score += 1.0;
                            } else {
                                switch (ftn_type) {
                                    case EXACT_FITNESS:
                                        break;
                                    case ATANG_FITNESS:
                                        part_scores.push(1.0 - 2.0 * atan((delta_err - m_attr_size) * m_ftn_scale) / PI);
                                        break;
                                    case INVER_FITNESS:
                                        part_scores.push(1.0 / (1.0 + (delta_err - m_attr_size) * m_ftn_scale));
                                        break;
                                    default:
                                        (void) throwException(env, IllegalArgumentException,
                                                "Unsupported fitness type!");
                                }
                            }
                        }
                    }

                    /**
                     * Allows to compute the error of the given controller.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     * @param scale scale for the case scaling is requested
                     * @param shift shift for the case scaling is requested
                     */
                    template<bool is_complex, fitness_type ftn_type>
                    void compute_fitness_scale(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        //Allocate the inputs array
                        double * inputs = new double[m_num_states];
                        //Set the cursor
                        double const * cursor = &m_all_data[0];
                        //Iterate over controller states and compute inputs
                        abs_type state_idx = 0;
                        double max_input = 0.0;
                        double min_input = DBL_MAX;
                        while (cursor < m_cursor_max) {
                            //Compute the candidate's input for the state
                            const double ctr_value = wrap.compute_input(cursor, m_ss_dim);
                            //Check if the value is sensible
                            if (!isnan(ctr_value) && !isinf(ctr_value)) {
                                //Store the value
                                inputs[state_idx] = ctr_value;
                                //Compute the min/max values
                                max_input = max(max_input, ctr_value);
                                min_input = min(min_input, ctr_value);
                                //Shift the cursor to where the inputs count is
                                cursor += m_ss_dim;
                                //Get the inputs count
                                const abs_type ips_cnt = static_cast<abs_type> (*cursor);
                                //Shift through the inputs part
                                cursor += 1 + m_is_dim * ips_cnt;
                                //Move to the next state index
                                ++state_idx;
                            } else {
                                //The value is off so it makes no sense to proceed
                                *this << "ERROR: Dof " << wrap.get_dof_idx()
                                        << " individual: " << wrap.get_name()
                                        << " input value is: " << ctr_value;
                                delete[] inputs;
                                return;
                            }
                        }

                        //Get the dof index relative to the largest state dof index
                        const int is_dof_idx = wrap.get_dof_idx();
                        //Initialize the priority queue for safe computations
                        scores_queue part_scores;
                        //Get the original controller min/max values
                        const pair<double, double> & min_max = m_min_max[is_dof_idx];
                        //Compute the scaling factor
                        const double orig_span = (min_max.second - min_max.first);
                        const double func_span = (max_input - min_input);

                        //Set the scaling and shifting values                        
                        scale = (func_span == 0.0) ? 1.0 : orig_span / func_span;
                        shift = min_max.first - min_input * scale;

                        //Re-initialize the cursor
                        cursor = &m_all_data[0];
                        //Initialize the scores
                        double part_score_sum = 0.0;
                        double exact_score = 0.0;
                        double part_score = 0.0;
                        //Iterate over actual controller and compute errors
                        state_idx = 0;
                        while (cursor < m_cursor_max) {
                            //Compute the state error
                            double delta_err = DBL_MAX;
                            const double adj_input = inputs[state_idx] * scale + shift;
                            if (compute_state_error(adj_input, is_dof_idx, cursor, delta_err)) {
                                //Compute the scores
                                update_scores<is_complex, ftn_type>(env, delta_err,
                                        exact_score, part_score, part_scores);
                            } else {
                                //We still got an error that is off, so we give up 
                                delete[] inputs;
                                return;
                            }
                            ++state_idx;
                        }

                        //Compute the resulting fitness values
                        while (!part_scores.empty()) {
                            part_score_sum += part_scores.top();
                            part_scores.pop();
                        }

                        exact_ftn = exact_score / m_num_states;
                        if (is_complex) {
                            const double ext_ftn = (part_score + part_score_sum) / m_num_states;
                            req_ftn = sqrt(exact_ftn * exact_ftn + ext_ftn * ext_ftn);
                        } else {
                            req_ftn = exact_ftn;
                        }
                        //Clear memory
                        delete[] inputs;
                    }

                    /**
                     * Allows to compute the error of the given controller.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     */
                    template<bool is_complex, fitness_type ftn_type>
                    void compute_fitness_plain(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn) {
                        //Get the dof index relative to the largest state dof index
                        const int is_dof_idx = wrap.get_dof_idx();
                        //Re-set the values to zero first
                        double exact_score = 0.0;
                        double part_score = 0.0;
                        //Create an inverted priority queue
                        scores_queue part_scores;

                        //Define the cursor
                        double const * cursor = &m_all_data[0];
                        double part_score_sum = 0.0;
                        while (cursor < m_cursor_max) {
                            //Compute the state error
                            double delta_err = DBL_MAX;
                            const double ctr_value = wrap.compute_input(cursor, m_ss_dim);
                            if (compute_state_error(ctr_value, is_dof_idx, cursor, delta_err)) {
                                //Compute the scores
                                update_scores<is_complex, ftn_type>(env, delta_err,
                                        exact_score, part_score, part_scores);
                            } else {
                                //We still got an error that is off, so we give up 
                                return;
                            }
                        }

                        //Compute the resulting fitness values
                        while (!part_scores.empty()) {
                            part_score_sum += part_scores.top();
                            part_scores.pop();
                        }

                        exact_ftn = exact_score / m_num_states;
                        if (is_complex) {
                            const double ext_ftn = (part_score + part_score_sum) / m_num_states;
                            req_ftn = sqrt(exact_ftn * exact_ftn + ext_ftn * ext_ftn);
                        } else {
                            req_ftn = exact_ftn;
                        }
                    }

                    /**
                     * Allows to store the unfit points into a BDD.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     */
                    void add_unfit_points_to_bdd(JNIEnv * env, ctrl_wrapper & wrap) {
                        //Get the dof index relative to the largest state dof index
                        const int is_dof_idx = wrap.get_dof_idx();
                        //The variables to store the number of fit and unfit points
                        long total = 0.0, unfit = 0.0;

                        //Define the cursor
                        double const * cursor = &m_all_data[0];
                        abs_type state_abs[m_ss_dim];
                        while (cursor < m_cursor_max) {
                            //Compute the state error
                            double delta_err = DBL_MAX;
                            const double ctr_value = wrap.compute_input(cursor, m_ss_dim);
                            double const * state_dbl = cursor;
                            //If there are nun of inf values or the error delta is
                            //larger than the maximum then this point is not fit.
                            if (!compute_state_error(ctr_value, is_dof_idx, cursor, delta_err)
                                    || (delta_err >= MAX_ATTRACTOR_SIZE)) {
                                //Copy the values into the state array
                                for (int idx = 0; idx < m_ss_dim; ++idx) {
                                    state_abs[idx] = *state_dbl;
                                    ++state_dbl;
                                }

                                //Convert the abstract state to BDD and add to the total
                                m_unfit_bdd = m_unfit_bdd | m_p_ss_mgr->x_to_bdd(state_abs);

                                //Increment the number of unfit points
                                ++unfit;
                            }
                            //Increment the number of fit points
                            ++total;
                        }

                        *this << "The number of unfit points for dof "
                        << wrap.get_dof_idx() << " is " << unfit << "/" << total;
                    }

                    /**
                     * Allows to convert the points array into the internal data structures
                     * @param num_pg the number of points in the array
                     * @param arr_gp the points array
                     * @return the number of states
                     */
                    abs_type convert_points_to_data() {
                        //Initialize the inputs and states managers
                        initialize_managers();

                        //Initialize the min-max pairs
                        m_min_max.clear();
                        for (abs_type idx = 0; idx < m_is_dim; ++idx) {
                            m_min_max.emplace_back(DBL_MAX, 0.0);
                        }

                        //Get the number the states with inputs
                        raw_data all_states = m_p_ss_mgr->get_points();
                        const int num_states = all_states.size() / m_ss_dim;
                        raw_data state(m_ss_dim), input(m_is_dim);
                        abs_data astate(m_ss_dim), ainput(m_is_dim);

                        //Iterate over the states, get the corresponding
                        //inputs and add them to the estimator set by ids
                        auto state_begin = all_states.begin();
                        while (state_begin != all_states.end()) {
                            //Get a new state vector
                            state.assign(state_begin, state_begin + m_ss_dim);

                            //Get the list of inputs
                            raw_data state_inputs = m_p_ctr->restriction(*m_p_cudd, *m_p_ctrl_bdd, state);

                            //Convert state to abstract state and add to data
                            m_p_ss_mgr->xtois(state, astate);
                            m_all_data.insert(m_all_data.end(), astate.begin(), astate.end());

                            //Get the number of inputs and add to the data
                            const int ips_cnt = state_inputs.size() / m_is_dim;
                            m_all_data.push_back(ips_cnt);

                            //Iterate over the inputs
                            auto input_begin = state_inputs.begin();
                            while (input_begin != state_inputs.end()) {
                                //Get a new input vector
                                input.assign(input_begin, input_begin + m_is_dim);

                                //Convert state to abstract state and add to data
                                m_p_is_mgr->xtois(input, ainput);
                                m_all_data.insert(m_all_data.end(), ainput.begin(), ainput.end());

                                //Compute the min and max values per input dimension
                                for (abs_type idx = 0; idx < m_is_dim; ++idx) {
                                    pair<double, double> & elem = m_min_max[idx];
                                    elem.first = min(elem.first, static_cast<double> (ainput[idx]));
                                    elem.second = max(elem.second, static_cast<double> (ainput[idx]));
                                }

                                //Move forward in the list of inputs
                                input_begin += m_is_dim;
                            }

                            //Move forward in the list of states
                            state_begin += m_ss_dim;
                        }

                        //Add the MAX_STATE_DEVIATION to the bounds
                        for (abs_type idx = 0; idx < m_is_dim; ++idx) {
                            m_min_max[idx].first -= MAX_STATE_DEVIATION;
                            m_min_max[idx].second += MAX_STATE_DEVIATION;
                        }

                        //Set the maximum cursor value
                        m_cursor_max = &m_all_data[m_all_data.size() - 1] + 1;

                        //Return the number of states
                        return num_states;
                    }

                    /**
                     * Allows to re-initialize the inputs and states managers
                     */
                    void initialize_managers() {
                        //Delete any previously present managers
                        if (m_p_ss_mgr) delete m_p_ss_mgr;
                        if (m_p_is_mgr) delete m_p_is_mgr;

                        //Initialize the input and state manager
                        m_p_is_mgr = new inputs_mgr(*m_p_ctr, m_ss_dim);
                        m_p_ss_mgr = new states_mgr(*m_p_ctr, m_ss_dim,
                                *m_p_ctrl_bdd, *m_p_cudd,
                                m_p_is_mgr->get_inputs_set());
                    }

                    /**
                     * Allows to set the state-space size
                     * NOT thread safe!
                     * @param env the JNI environment
                     * @param size the state-space size
                     */
                    void set_ss_dim(JNIEnv * const env, const jint ss_dim) {
                        if (ss_dim < 1 || ss_dim >= m_p_ctr->get_dim()) {
                            (void) throwException(env, IllegalArgumentException,
                                    "The state-space size must be [1, #dofs)");
                        } else {
                            //if the states were not computed yet or the state 
                            //space size has changed then we need to re-compute
                            if (m_all_data.empty() || (m_ss_dim != ss_dim)) {
                                //Store the size
                                m_ss_dim = ss_dim;
                                m_is_dim = (m_p_ctr->get_dim() - m_ss_dim);
                                //Clear the old data
                                m_all_data.clear();
                                *this << "Start extracting grid points";
                                //Convert the points into the internal data structures
                                m_num_states = convert_points_to_data();
                                *this << "The state-space size is = " << m_num_states;
                            }
                        }
                    }

                    /**
                     * Allows to set the attractor size, must be <= 0.5
                     * @param env the JNI environment
                     * @param attr_size the attractor size
                     */
                    void set_attr_size(JNIEnv * env, const jdouble attr_size) {
                        if (m_is_complex) {
                            if (attr_size < 0 || attr_size >= MAX_ATTRACTOR_SIZE) {
                                (void) throwException(env, IllegalArgumentException,
                                        "Improper attractor value, must be within [0, 0.5)!");
                            } else {
                                m_attr_size = attr_size;
                            }
                        }
                    }

                    /**
                     * Allows to get the Java object field value
                     * @param env the JNI environment
                     * @param cls the object class
                     * @param obj the object itself
                     * @param fid_name the field name
                     * @return the field value
                     */
                    int get_int_field(JNIEnv * env, jclass cls, jobject obj,
                            const char * fid_name) {
                        jfieldID fid = env->GetFieldID(cls, fid_name, "I");
                        return env->GetIntField(obj, fid);
                    }

                    /**
                     * Allows to get the Java object field value
                     * @param env the JNI environment
                     * @param cls the object class
                     * @param obj the object itself
                     * @param fid_name the field name
                     * @return the field value
                     */
                    double get_double_field(JNIEnv * env, jclass cls, jobject obj,
                            const char * fid_name) {
                        jfieldID fid = env->GetFieldID(cls, fid_name, "D");
                        return env->GetDoubleField(obj, fid);
                    }

                    /**
                     * Allows to get the Java object field value
                     * @param env the JNI environment
                     * @param cls the object class
                     * @param obj the object itself
                     * @param fid_name the field name
                     * @return the field value
                     */
                    bool get_bool_field(JNIEnv * env, jclass cls, jobject obj,
                            const char * fid_name) {
                        jfieldID fid = env->GetFieldID(cls, fid_name, "Z");
                        return env->GetBooleanField(obj, fid);
                    }

                    /**
                     * Allows to compute the error of the given controller.
                     * @param env the JNI environment
                     * @param wrap the controller's wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     * @param scale scale for the case scaling is requested
                     * @param shift shift for the case scaling is requested
                     */
                    template<bool IS_COMPLEX>
                    void compute_fitness(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        switch (m_ftn_type) {
                            case EXACT_FITNESS:
                                if (m_is_scale) {
                                    compute_fitness_scale< IS_COMPLEX, EXACT_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn, scale, shift);
                                } else {
                                    compute_fitness_plain< IS_COMPLEX, EXACT_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn);
                                }
                                break;
                            case ATANG_FITNESS:
                                if (m_is_scale) {
                                    compute_fitness_scale< IS_COMPLEX, ATANG_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn, scale, shift);
                                } else {
                                    compute_fitness_plain< IS_COMPLEX, ATANG_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn);
                                }
                                break;
                            case INVER_FITNESS:
                                if (m_is_scale) {
                                    compute_fitness_scale< IS_COMPLEX, INVER_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn, scale, shift);
                                } else {
                                    compute_fitness_plain< IS_COMPLEX, INVER_FITNESS >
                                            (env, wrap, exact_ftn, req_ftn);
                                }
                                break;
                            default:
                                (void) throwException(env, IllegalArgumentException,
                                        "Unsupported fitness type!");
                        }
                    }

                public:

                    ftn_computer() :
                    m_log_file(), m_p_cudd(NULL), m_p_ctrl_bdd(NULL),
                    m_p_ctr(NULL), m_unfit_bdd(), m_p_is_mgr(NULL),
                    m_p_ss_mgr(NULL), m_ss_dim(0), m_is_dim(0),
                    m_num_states(0), m_all_data(), m_cursor_max(NULL),
                    m_attr_size(0.0), m_min_max(), m_is_scale(false),
                    m_ftn_type(UNDEF_FITNESS), m_is_complex(true),
                    m_ftn_scale(1.0) {
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
                        string log_file_name = string(file_name) + ".sr.log";

                        //(Re-)Open the log file
                        if (m_log_file.is_open()) {
                            m_log_file.close();
                        }
                        m_log_file.open(log_file_name);
                        m_log_file.clear();

                        *this << "Opened the log-file: " << log_file_name;

                        //Read controller from file
                        *this << "Start initializing data structures";

                        //De-couple the BDD from the CUDD manager
                        m_unfit_bdd = BDD();

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
                        *this << "Start loading the controller: " << file_name;

                        //Load the controller
                        bool is_fail = false;
                        if (!read_from_file(*m_p_cudd, *m_p_ctr, *m_p_ctrl_bdd, file_name)) {
                            is_fail = true;
                        }
                        *this << "Loading: " << file_name << " is finished";

                        //Throw if could not load otherwise extract the bdd data
                        if (is_fail) {
                            (void) throwException(env, FileNotFoundException,
                                    "Failed to read the controller file!");
                        }

                        return m_p_ctr->get_dim();
                    }

                    /**
                     * Allows to get the number of points of the state-space grid.
                     * @return the number of points of the state-space grid
                     */
                    int get_state_space_size(JNIEnv * env) {
                        *this << "Retrieving the number of points on the state-space grid";
                        int size = 1;
                        if (m_p_ctr != NULL) {
                            const vector<abs_type> nu_gp = m_p_ctr->get_no_gp_per_dim();
                            for (int idx = 1; idx < m_p_ctr->get_dim(); ++idx) {
                                size *= nu_gp[idx];
                            }
                        } else {
                            (void) throwException(env, IllegalStateException,
                                    "The controller is not loaded yet!");
                        }
                        return size;
                    }

                    /**
                     * Allows to configure the fitness computer based on the configuration object
                     * @param env the JNI environment
                     * @param obj the configuration object
                     */
                    void configure(JNIEnv * env, jobject obj) {
                        jclass cls = env->GetObjectClass(obj);
                        const int ss_size = get_int_field(env, cls, obj, "m_ss_size");
                        set_ss_dim(env, ss_size);
                        m_ftn_type = (fitness_type) get_int_field(env, cls, obj, "m_ftn_type");
                        const double attr_size = get_double_field(env, cls, obj, "m_attr_size");
                        m_is_scale = get_bool_field(env, cls, obj, "m_is_scale");
                        m_is_complex = get_bool_field(env, cls, obj, "m_is_complex");
                        m_ftn_scale = get_double_field(env, cls, obj, "m_ftn_scale");
                        if (m_ftn_scale == 0.0) {
                            m_ftn_scale = 1.0;
                        }
                        set_attr_size(env, attr_size);
                    }

                    /**
                     * Allows to start the unfit points export by re-setting
                     * the internal data storage for the unfit points
                     * @param env the environment
                     */
                    void start_unfit_export(JNIEnv * env) {
                        *this << "Starting unfit points export";
                        m_unfit_bdd = m_p_cudd->bddZero();
                    }

                    /**
                     * Allows to add the unfit points of this individual into
                     * the unfit points to be exported into the BDD
                     * @param env the JNI environment
                     * @param wrap the controller wrapper
                     */
                    void compute_unfit_points(JNIEnv * env, ctrl_wrapper & wrap) {
                        *this << "Starting unfit points export of " << wrap.get_name()
                                << " for dof " << wrap.get_dof_idx();
                        if (m_p_ctr != NULL) {
                            if (m_ss_dim > 0) {
                                const int dof_idx = wrap.get_dof_idx();
                                if ((dof_idx >= 0) && (dof_idx < (m_p_ctr->get_dim() - m_ss_dim))) {
                                    add_unfit_points_to_bdd(env, wrap);
                                } else {
                                    (void) throwException(env, IllegalArgumentException,
                                            "Improper input space dimension index!");
                                }
                            } else {
                                (void) throwException(env, IllegalStateException,
                                        "The state-space size is not set!");
                            }
                        } else {
                            (void) throwException(env, IllegalStateException,
                                    "The controller is not loaded yet!");
                        }
                    }

                    /**
                     * Allows to store the collected so-far unfit points into the BDD.
                     * @param env the JNI environment
                     * @param file_name the BDD file name
                     */
                    void finish_unfit_export(JNIEnv * env, const char * file_name) {
                        *this << "Dumping unfit points into " << file_name;
                        store_controller(*m_p_cudd,
                                m_p_ss_mgr->get_states_set(),
                                m_unfit_bdd, string(file_name));
                    }

                    /**
                     * Allows to compute fitness of the given controller wrapper
                     * @param env the JNI environment
                     * @param wrap the controller wrapper
                     * @param exact_ftn the exact fitness
                     * @param req_ftn the requested fitness
                     * @param scale scale for the case scaling is requested
                     * @param shift shift for the case scaling is requested
                     */
                    void compute(JNIEnv * env, ctrl_wrapper & wrap,
                            double & exact_ftn, double & req_ftn,
                            double & scale, double & shift) {
                        if (m_p_ctr != NULL) {
                            if (m_ss_dim > 0) {
                                const int dof_idx = wrap.get_dof_idx();
                                if ((dof_idx >= 0) && (dof_idx < (m_p_ctr->get_dim() - m_ss_dim))) {
                                    if (m_is_complex) {
                                        compute_fitness<true>(env, wrap, exact_ftn, req_ftn, scale, shift);
                                    } else {
                                        compute_fitness<false>(env, wrap, exact_ftn, req_ftn, scale, shift);
                                    }
                                } else {
                                    (void) throwException(env, IllegalArgumentException,
                                            "Improper input space dimension index!");
                                }
                            } else {
                                (void) throwException(env, IllegalStateException,
                                        "The state-space size is not set!");
                            }
                        } else {
                            (void) throwException(env, IllegalStateException,
                                    "The controller is not loaded yet!");
                        }
                    }

                    /**
                     * The basic destructor
                     */
                    virtual ~ftn_computer() {
                        if (m_log_file.is_open()) {
                            m_log_file.close();
                        }
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
                        //De-couple the BDD from the cudd manager
                        m_unfit_bdd = BDD();
                        if (m_p_cudd) {
                            delete m_p_cudd;
                        }
                    }
                };

                ostream& operator<<(ftn_computer& strm, const char * str) {
                    time_t now = time(0);
                    string time(ctime(&now));
                    time = time.substr(0, time.length() - 1);
                    return strm.m_log_file << time << "\b: " << str
                            << std::endl << std::flush;
                }
            }
        }
    }
}

#endif /* FTN_COMPUTER_HH */

