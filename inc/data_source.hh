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

#include <mutex>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#define SCOTS_BDD 1
#include "scots.hh"

#include "inputs_mgr.hh"
#include "states_mgr.hh"
#include "input_output.hh"
#include "jni_throw.hh"

#include "info_logger.hh"
#include "config_wrapper.hh"
#include "random_hypercube.hh"
#include "var_bisector.hh"

using namespace std;
using namespace scots;
using namespace tud::ctrl::scots::optimal;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                typedef vector<double> raw_data;
                typedef vector<abs_type> abs_data;
                typedef unordered_set<uint64_t> sample_set;

                /**
                 * This class is used to load and store the controller's data
                 */
                class data_source {
                private:
                    //Stores the mutex for thread safety
                    mutable mutex m_mutex;
                    //Will store the state
                    mutable raw_data m_state;
                    //Stores the cached inputs
                    mutable unordered_map<abs_type, raw_data> m_cache;

                    //Stores the number of the domain elements
                    int64_t m_dom_size;

                    //Stores the sample to use for fitting
                    raw_data m_sample_data;
                    //Stores the data sample size
                    abs_type m_sample_size;

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

                    //Stores the min/max values for the input dofs
                    vector<pair<double, double>> m_min_max;

                    /**
                     * The basic constructor.
                     */
                    data_source()
                    : m_state(), m_cache(), m_p_cudd(NULL), m_p_ctrl_bdd(NULL),
                    m_p_ctr(NULL), m_p_is_mgr(NULL), m_p_ss_mgr(NULL),
                    m_dom_size(0), m_sample_data(), m_sample_size(), m_min_max() {
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
                        LOG("Start initializing data structures");

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
                        LOG("Start loading the controller: " << file_name);

                        //Load the controller
                        bool is_fail = false;
                        if (!read_from_file(*m_p_cudd, *m_p_ctr, *m_p_ctrl_bdd, file_name)) {
                            is_fail = true;
                        }
                        LOG("Loading: " << file_name << " is finished");

                        //Throw if could not load otherwise extract the bdd data
                        if (is_fail) {
                            (void) throwException(env, FileNotFoundException,
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
                        LOG("Retrieving the number of points on the state-space grid");
                        int size = 1;
                        if (m_p_ctr != NULL) {
                            const vector<abs_type> nu_gp = m_p_ctr->get_no_gp_per_dim();
                            for (int idx = 0; idx < ss_dim; ++idx) {
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
                     * @param cfg the configuration object
                     */
                    void configure(JNIEnv * env, jobject cfg) {
                        //Create a new configuration wrapper
                        config_wrapper config(env, cfg, m_p_ctr->get_dim());

                        //Check that the configuration is correct
                        config.finalize(env);

                        //Store the reference to the configuration object
                        m_config = config;

                        //Clean the old data if any
                        re_initialize_data();

                        //Convert the points into the internal data structures, only if needed!
                        //Note the the Monte Carlo fitness (MC and RSS) do not support scaling.
                        if (m_config.m_is_mc) {
                            if (m_config.m_is_scale) {
                                prepare_sim_data<true>();
                            } else {
                                prepare_sim_data<false>();
                            }
                        } else {
                            if (m_config.m_is_scale) {
                                prepare_num_data<true>();
                            } else {
                                prepare_num_data<false>();
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
                        //If we did the Monte-Carlo simulations then the
                        //plain data is not complete clear it and re-fill
                        //For numeric computations the plain data must
                        //have been fully loaded.
                        if (m_config.m_is_mc) {
                            //Clear the sample data
                            clear_sample_data();
                            //Re-fill the sample with complete data,
                            //no min/max computations are needed
                            prepare_num_data<false>();
                        }
                    }

                    /**
                     * Allows to compute the restriction of the state in a thread-safe way.
                     * @param state the state to restrict to
                     * @return the vector of resulting inputs
                     */
                    inline const raw_data & restriction(const uint64_t id) const {
                        lock_guard<mutex> guard(m_mutex);

                        //Check if the state is cached
                        auto val_iter = m_cache.find(id);
                        if (val_iter != m_cache.end()) {
                            return val_iter->second;
                        } else {
                            //Compute the state
                            m_p_ss_mgr->itox(id, m_state);

                            //Compute the restriction                        
                            auto result = m_cache.emplace(id,
                                    m_p_ctr->restriction(*m_p_cudd, *m_p_ctrl_bdd, m_state));

                            return result.first->second;
                        }
                    }

                    /**
                     * Allows to retrieve the sample size
                     * @return the sample size
                     */
                    inline abs_type get_sample_size() const {
                        return m_sample_size;
                    }

                    /**
                     * Allows to get the start of the plain data
                     * @return the start of the plain data
                     */
                    inline double const * get_sample_start() const {
                        return &m_sample_data[0];
                    }

                    /**
                     * Allows to get the end of the plain data
                     * @return the end of the plain data
                     */
                    inline double const * get_sample_end() const {
                        return &m_sample_data[m_sample_data.size() - 1] + 1;
                    }

                private:
                    //Stores the configuration object
                    config_wrapper m_config;

                    /**
                     * Allows to clear the sample data
                     * @param is_final if true then this is a final clean-up
                     */
                    void clear_sample_data() {
                        m_sample_data.clear();
                        m_sample_size = 0;
                    }

                    /**
                     * Clears the old data and does initial re-initialization
                     */
                    void re_initialize_data() {
                        //Clear the cached data
                        m_cache.clear();

                        //Initialize the min-max pairs
                        m_min_max.clear();

                        //Delete any previously present managers
                        if (m_p_ss_mgr) delete m_p_ss_mgr;
                        if (m_p_is_mgr) delete m_p_is_mgr;

                        //Initialize the input and state manager
                        m_p_is_mgr = new inputs_mgr(*m_p_ctr, m_config.m_num_ss_dim);
                        m_p_ss_mgr = new states_mgr(*m_p_ctr, m_config.m_num_ss_dim,
                                *m_p_ctrl_bdd, *m_p_cudd,
                                m_p_is_mgr->get_inputs_set());

                        //Clear the old data points
                        clear_sample_data();

                        //Set the number of states value
                        m_dom_size = m_p_ss_mgr->get_size();
                        LOG("The controller domain size: " << m_dom_size);
                    }

                    /**
                     * Pre-initialize the min max values if we need scaling
                     */
                    template<bool IS_SCALE_FIT>
                    void pre_set_min_max() {
                        //Initialize the min/max values
                        if (IS_SCALE_FIT) {
                            for (abs_type idx = 0; idx < m_config.m_num_is_dim; ++idx) {
                                m_min_max.emplace_back(DBL_MAX, 0.0);
                            }
                        }
                    }

                    /**
                     * Allows to update min/max input values if scaling is needed
                     * @param ainput the input vector
                     */
                    template<bool IS_SCALE_FIT>
                    void update_min_max_values(const abs_data & ainput) {
                        //Compute the min and max values per input dimension
                        if (IS_SCALE_FIT) {
                            for (abs_type idx = 0; idx < m_config.m_num_is_dim; ++idx) {
                                pair<double, double> & elem = m_min_max[idx];
                                elem.first = min(elem.first, static_cast<double> (ainput[idx]));
                                elem.second = max(elem.second, static_cast<double> (ainput[idx]));
                            }
                        }
                    }

                    /**
                     * Post-process the min max values if we need scaling
                     */
                    template<bool IS_SCALE_FIT>
                    void post_process_min_max() {
                        //Add the MAX_STATE_DEVIATION to the bounds
                        if (IS_SCALE_FIT) {
                            for (abs_type idx = 0; idx < m_config.m_num_is_dim; ++idx) {
                                m_min_max[idx].first -= MAX_STATE_DEVIATION;
                                m_min_max[idx].second += MAX_STATE_DEVIATION;
                            }
                        }
                    }

                    /**
                     * Allows to store the valid domain state data
                     * @param state_abs the abstract state vector representation
                     * @param state_dbl the abstract state vector representation as a vector of doubles
                     * @param state the container for storing the state
                     * @param inputs the temporary container for the state inputs
                     * @param call_back the call back object to use if the sample is added
                     * @param p_bisect the pointer to the bisector, if IS_BISECTOR is true or NULL, default to NULL
                     * @return true if this is a valid domain state and thus was stored
                     */
                    template<bool IS_SCALE_FIT, bool IS_BISECTOR>
                    bool store_dom_state_data(uint64_t const * state_abs,
                            double const * state_dbl, raw_data & state,
                            raw_data & inputs, var_bisector * p_bisect = NULL) {
                        //Compute the actual state
                        m_p_ss_mgr->Itox(state_abs, state);

                        //Compute the restriction                        
                        m_p_ctr->restriction(*m_p_cudd, *m_p_ctrl_bdd, state, inputs);

                        //Register the data if it from the domain, i.e. there are inputs
                        return register_plain_data<IS_SCALE_FIT, IS_BISECTOR>(
                                state, inputs, p_bisect);
                    }

                    /**
                     * Allows to check if this is a valid domain state and cache it and its inputs
                     * @param state_abs the abstract state vector representation
                     * @param state_dbl the abstract state vector representation as a vector of doubles
                     * @param samples the set storing all the unique sample element ids
                     * @param state the container for storing the state
                     * @param inputs the temporary container for the state inputs
                     * @param p_bisect the pointer to the bisector, if IS_BISECTOR is true or NULL, default to NULL
                     * @return true if this is a new valid domain state, otherwise false
                     */
                    template<bool IS_SCALE_FIT, bool IS_BISECTOR>
                    bool cache_sample_element(uint64_t const * state_abs,
                            double const * state_dbl, sample_set & samples,
                            raw_data & state, raw_data & inputs,
                            var_bisector * p_bisect = NULL) {
                        bool is_valid = false;
                        abs_type state_id;
                        //Check if the point is on the grid (this is more a safety check)
                        if (m_p_ss_mgr->istoi(state_abs, state_id)) {
                            //Check if the point is not cached yet then try it out
                            auto res = samples.insert(state_id);
                            if (res.second) {
                                //Store the point if it is in the domain
                                is_valid = store_dom_state_data<IS_SCALE_FIT, IS_BISECTOR>(
                                        state_abs, state_dbl, state, inputs, p_bisect);
                                //Check if the point turned out to be 
                                //not from the domain then remove it
                                if (!is_valid) {
                                    samples.erase(res.first);
                                }
                            }
                        }
                        return is_valid;
                    }

                    /**
                     * Allows to prepare data for the Monte-Carlo fitness
                     * @param cube the random hypercube to sample from
                     * @param sample_size the desired sample size
                     * @param samples the set to store sample state ids
                     * @param state_abs internal temporary container for state abstract vector
                     * @param state_dbl internal temporary container for state abstract vector as double
                     * @param state internal temporary container for real state vector
                     * @param inputs internal temporary container for real input vectors
                     * @param p_bisect the pointer to the bisector, if IS_BISECTOR is true or NULL, default to NULL
                     */
                    template<bool IS_SCALE_FIT, bool IS_BISECTOR>
                    void prepare_sim_mc_data(
                            const int64_t sample_size,
                            random_hypercube & cube, sample_set & samples,
                            uint64_t * state_abs, double * state_dbl,
                            raw_data & state, raw_data & inputs,
                            var_bisector * p_bisect = NULL) {
                        LOG("Monte-Carlo simulations, sample size: " << sample_size << " "
                                << (IS_BISECTOR ? "WITH" : "WITHOUT") << " bisector");

                        //Prepare samples
                        size_t ext_smp_cnt = 0;
                        while (ext_smp_cnt < sample_size) {
                            //Generate a random abstract state
                            cube.get_random(state_abs, state_dbl);

                            //Cache the domain state with its inputs
                            if (cache_sample_element<IS_SCALE_FIT, IS_BISECTOR>(
                                    state_abs, state_dbl, samples, state,
                                    inputs, p_bisect)) {
                                //Increment the sample count
                                ++ext_smp_cnt;
                            }
                        }
                    }

                    /**
                     * Allows to bisect the current hypercube and to Recursion into the sub-cubes
                     * @param sample_size the desired sample size
                     * @param bisect the bisector with some pre-sampled data for bisection
                     * @param cube the random hypercube to sample from
                     * @param samples the set to store sample state ids
                     * @param state_abs internal temporary container for state abstract vector
                     * @param state_dbl internal temporary container for state abstract vector as double
                     * @param state internal temporary container for real state vector
                     * @param inputs internal temporary container for real input vectors
                     */
                    template<bool IS_SCALE_FIT>
                    void prepare_sim_rss_data(
                            const int64_t sample_size,
                            var_bisector & bisect, random_hypercube & cube,
                            sample_set & samples, uint64_t * state_abs,
                            double * state_dbl, raw_data & state,
                            raw_data & inputs) {
                        LOG("Remaining sample size: " << sample_size);

                        //If the bisection dof could be found then split
                        const int dim_idx = bisect.get_bisect_dof();
                        LOG("The bisection dof idx: " << dim_idx);
                        if ((dim_idx >= 0) && (dim_idx < m_config.m_num_ss_dim)) {
                            //Get the bisection ratio
                            const double ratio = bisect.get_lr_ratio();

                            LOG("The bisection ratio is: " << ratio);

                            //Get the sub hypercubes as bisected
                            random_hypercube * p_lc = NULL, * p_rc = NULL;
                            uint64_t lr = 0, rl = 0;
                            bisect.get_lr_rl_dounds(lr, rl);
                            cube.split(dim_idx, lr, rl, p_lc, p_rc);

                            //Recursively go to left hypercube
                            const int64_t lcs_size = sample_size * ratio;
                            prepare_sim_rss_data<IS_SCALE_FIT>(
                                    lcs_size, *p_lc, samples, state_abs,
                                    state_dbl, state, inputs);
                            delete p_lc;

                            //Recursively go to right hypercube
                            const int64_t rcs_size = sample_size - lcs_size;
                            prepare_sim_rss_data<IS_SCALE_FIT>(
                                    rcs_size, *p_rc, samples, state_abs,
                                    state_dbl, state, inputs);
                            delete p_rc;
                        } else {
                            //The bisection is not needed, just proceed with regular MC
                            prepare_sim_mc_data<IS_SCALE_FIT, false>(
                                    sample_size, cube, samples, state_abs,
                                    state_dbl, state, inputs);
                        }
                    }

                    /**
                     * Allows to prepare data for the Monte-Carlo Recursive Stratified Sampling fitness
                     * @param sample_size the desired sample size
                     * @param cube the random hypercube to sample from
                     * @param samples the set to store sample state ids
                     * @param state_abs internal temporary container for state abstract vector
                     * @param state_dbl internal temporary container for state abstract vector as double
                     * @param state internal temporary container for real state vector
                     * @param inputs internal temporary container for real input vectors
                     */
                    template<bool IS_SCALE_FIT>
                    void prepare_sim_rss_data(
                            const int64_t sample_size,
                            random_hypercube & cube, sample_set & samples,
                            uint64_t * state_abs, double * state_dbl,
                            raw_data & state, raw_data & inputs) {
                        LOG("Remaining sample size: " << sample_size);

                        //Compute the bisection sample size
                        const int64_t bs_size = m_config.m_bisect_ratio*sample_size;
                        //Check if we can proceed bisecting
                        if (bs_size >= m_config.m_bisect_size) {
                            //Declare the bisector
                            var_bisector bisect(cube, m_p_is_mgr->get_dim());

                            //Do the simple MC generation with the new sample
                            //registration inside the bisector
                            prepare_sim_mc_data<IS_SCALE_FIT, true>(
                                    bs_size, cube, samples, state_abs,
                                    state_dbl, state, inputs, &bisect);

                            //Indicate that the sampling is finished
                            bisect.done_sampling();

                            //Bisect the hypercube based on the bisector statistics
                            //Notice that the sample size is to be reduced by the
                            //number of bisection samples used to decide on bisection
                            prepare_sim_rss_data<IS_SCALE_FIT>(
                                    sample_size - bs_size, bisect, cube, samples,
                                    state_abs, state_dbl, state, inputs);
                        } else {
                            //The bisection is not needed, just proceed with regular MC
                            prepare_sim_mc_data<IS_SCALE_FIT, false>(
                                    sample_size, cube, samples, state_abs,
                                    state_dbl, state, inputs);
                        }
                    }

                    /**
                     * Allows to prepare data for the Monte-Carlo fitness
                     */
                    template<bool IS_SCALE_FIT>
                    void prepare_sim_data() {
                        //Update the sample size, it can not be larger
                        //than the number of domain elements
                        const int64_t sample_size = min(m_config.m_sample_size, m_dom_size);
                        LOG("Requested sample size: " << m_config.m_sample_size <<
                                ", re-computed sample size: " << sample_size);

                        //Initialize the random hypercube
                        uint64_t * ll = new uint64_t[m_config.m_num_ss_dim];
                        uint64_t * ur = new uint64_t[m_config.m_num_ss_dim];
                        vector<abs_type> num_gp_per_dof = m_p_ss_mgr->get_no_gp_per_dim();
                        for (int idx = 0; idx < m_config.m_num_ss_dim; ++idx) {
                            ll[idx] = 0;
                            ur[idx] = num_gp_per_dof[idx] - 1;
                        }
                        random_hypercube cube(m_config.m_num_ss_dim, ll, ur);

                        //Pre-set min/max values
                        pre_set_min_max<IS_SCALE_FIT>();

                        //Declare data containers
                        uint64_t state_abs[m_config.m_num_ss_dim];
                        double state_dbl[m_config.m_num_ss_dim];
                        raw_data state(m_config.m_num_ss_dim), inputs;

                        //Compute sample data for plain MC or RSS
                        sample_set samples;
                        if (m_config.m_is_rss) {
                            prepare_sim_rss_data<IS_SCALE_FIT>(
                                    sample_size, cube, samples, state_abs,
                                    state_dbl, state, inputs);
                        } else {
                            prepare_sim_mc_data<IS_SCALE_FIT, false>(
                                    sample_size, cube, samples, state_abs,
                                    state_dbl, state, inputs);
                        }
                        m_sample_size = samples.size();

                        LOG("Re-computed sample size: " << sample_size
                                << ", actual sample size: " << m_sample_size);

                        //Post-process min/max values
                        post_process_min_max<IS_SCALE_FIT>();
                    }

                    /**
                     * Allows to register the state inputs data
                     * @param state the actual state
                     * @param inputs the list of the state's actual inputs
                     * @param p_bisect the pointer to the bisector, if IS_BISECTOR is true or NULL, default to NULL
                     * @return true if the state has inputs
                     */
                    template<bool IS_SCALE_FIT, bool IS_BISECTOR>
                    bool register_plain_data(
                            const raw_data & state,
                            const raw_data & inputs,
                            var_bisector * p_bisect = NULL) {
                        //Define the result
                        bool result = true;

                        //Define state and abstract state containers
                        raw_data input(m_config.m_num_is_dim);
                        abs_data astate(m_config.m_num_ss_dim);
                        abs_data ainput(m_config.m_num_is_dim);

                        //Check if there are inputs
                        if (inputs.size() > 0) {
                            //If we need to store all data explicitly for numeric fitness computations
                            //Convert state to abstract state and add to data
                            m_p_ss_mgr->xtois(state, astate);
                            m_sample_data.insert(m_sample_data.end(), astate.begin(), astate.end());

                            //Indicate the beginning of a new sample
                            if (IS_BISECTOR) {
                                p_bisect->start_sample(astate);
                            }

                            //Get the number of inputs and add to the data
                            const int ips_cnt = inputs.size() / m_config.m_num_is_dim;
                            m_sample_data.push_back(ips_cnt);

                            //Iterate over the inputs
                            auto input_begin = inputs.begin();
                            while (input_begin != inputs.end()) {
                                //Get a new input vector
                                input.assign(input_begin, input_begin + m_config.m_num_is_dim);

                                //Convert state to abstract state and add to data
                                m_p_is_mgr->xtois(input, ainput);
                                //If we need to store all data explicitly for numeric fitness computations
                                m_sample_data.insert(m_sample_data.end(), ainput.begin(), ainput.end());

                                //Add a new state input to the bisector, if needed
                                if (IS_BISECTOR) {
                                    abs_type input_id = m_p_is_mgr->xtoi(ainput);
                                    p_bisect->add_input(input_id, ainput);
                                }

                                //Update the min/max values
                                update_min_max_values<IS_SCALE_FIT>(ainput);

                                //Move forward in the list of inputs
                                input_begin += m_config.m_num_is_dim;
                            }

                            //Indicate the end of a the sample
                            if (IS_BISECTOR) {
                                p_bisect->stop_sample();
                            }
                        } else {
                            result = false;
                        }

                        return result;
                    }

                    /**
                     * Allows to convert the points array into the internal data structures
                     */
                    template<bool IS_SCALE_FIT>
                    void prepare_num_data() {
                        LOG("Start extracting grid points");

                        //Pre-set min/max values
                        pre_set_min_max<IS_SCALE_FIT>();

                        //Get the number the states with inputs
                        raw_data states;
                        m_p_ss_mgr->get_points(states);

                        //Define state and abstract state containers
                        raw_data state(m_config.m_num_ss_dim), inputs;

                        LOG("Starting to extract " << m_dom_size << " domain states...");

                        //Iterate over the states, get the corresponding
                        //inputs and add them to the estimator set by ids
                        auto state_begin = states.begin();
                        int64_t new_cnt = 0, old_cnt = 0, state_cnt = 0;
                        while (state_begin != states.end()) {
                            //Get a new state vector
                            state.assign(state_begin, state_begin + m_config.m_num_ss_dim);

                            //Get the list of inputs
                            m_p_ctr->restriction(*m_p_cudd, *m_p_ctrl_bdd, state, inputs);

                            //Register data
                            register_plain_data<IS_SCALE_FIT, false>(state, inputs);

                            //Move forward in the list of states
                            state_begin += m_config.m_num_ss_dim;

                            //De the logging for convenience
                            new_cnt = state_cnt * 1000 / m_dom_size;
                            if (new_cnt != old_cnt) {
                                LOG("Extracted " << ((double) new_cnt) / 10.0
                                        << "% of the domain states.");
                                old_cnt = new_cnt;
                            }
                            ++state_cnt;
                        }

                        m_sample_size = m_p_ss_mgr->get_size();
                        LOG("The data sample size: " << m_sample_size);

                        //Post-process min/max values
                        post_process_min_max<IS_SCALE_FIT>();
                    }
                };
            }
        }
    }
}

#endif /* DATA_SOURCE_HH */

