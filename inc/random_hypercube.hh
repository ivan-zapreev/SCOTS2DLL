/* 
 * File:   random_hypercube.hh
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
 * Created on April 25, 2018, 3:18 PM
 */

#ifndef RANDOM_HYPERCUBE_HH
#define RANDOM_HYPERCUBE_HH

#include <mutex>
#include <cmath>
#include <cstring>
#include <random>
#include <gsl/gsl_qrng.h>

#include "jni_throw.hh"
#include "info_logger.hh"

using namespace std;

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                //Get the random generator to generate boolean values
                static auto bool_gen = bind(uniform_int_distribution<>(0, 1), default_random_engine());

                /**
                 * Represents a hypercube within which values can be generated 
                 * using pseudo random numbers. This class is thread safe!
                 */
                class random_hypercube {
                private:
                    //Stores the mutex for thread safety
                    mutex m_mutex;
                    //Stores the number of hypercube dimensions
                    const int m_num_dim;
                    //Stores the lower left coordinates
                    const uint64_t* m_ll;
                    //Stores the upper right coordinates
                    const uint64_t* m_ur;
                    //The array storing the raw point value, for internal use
                    double * m_raw_point;
                    //Stores an instance of a quasi random sequence
                    gsl_qrng * m_qr_seq;
                    //Flag indicates whether this instance owns the quasi random sequence generator
                    const bool m_is_qr_seq_own;
                    //Stores the scaling factors per dof
                    double* m_scale;
                    //Stores the shifting factors per dof
                    double* m_shift;

                public:

                    /**
                     * The basic constructor
                     * @param qr_seq the pointer to pre-allocated quasi random sequence
                     * @param is_qr_seq_own true if this object will own the quasi random 
                     * sequence and will be responsible for its deallocation
                     * @param num_dim the number of hypercube dimensions
                     * @param ll the lower left coordinates, not null
                     * @param ur the upper right coordinates, not null
                     */
                    random_hypercube(gsl_qrng * qr_seq, const bool is_qr_seq_own,
                            const int num_dim, const uint64_t* ll, const uint64_t* ur)
                    : m_mutex(), m_num_dim(num_dim), m_ll(ll), m_ur(ur),
                    m_qr_seq(qr_seq), m_is_qr_seq_own(is_qr_seq_own) {
                        //Allocate the raw point container
                        m_raw_point = new double[m_num_dim];
                        //Allocate the scaling and shifting factors
                        m_scale = new double[m_num_dim];
                        m_shift = new double[m_num_dim];
                        //Initialize the scaling and shifting factors
                        for (int dof_idx = 0; dof_idx < m_num_dim; ++dof_idx) {
                            //LOG("DOF[" << dof_idx << "] range [" << m_ll[dof_idx] << ", " << m_ur[dof_idx] << "]");
                            m_shift[dof_idx] = ((double) m_ll[dof_idx]) - 0.5;
                            m_scale[dof_idx] = (double) (m_ur[dof_idx] - m_ll[dof_idx] + 1);
                        }
                    }

                    /**
                     * The basic constructor
                     * @param seq_type the quasi random sequence type
                     * @param num_dim the number of hypercube dimensions
                     * @param ll the lower left coordinates, not null
                     * @param ur the upper right coordinates, not null
                     */
                    random_hypercube(const gsl_qrng_type * seq_type,
                            const int num_dim, const uint64_t* ll, const uint64_t* ur)
                    : random_hypercube(gsl_qrng_alloc(seq_type, num_dim),
                    true, num_dim, ll, ur) {
                    }

                    /**
                     * The basic constructor, uses Sobol Quasi Random Sequence
                     * @param num_dim the number of hypercube dimensions
                     * @param ll the lower left coordinates, not null
                     * @param ur the upper right coordinates, not null
                     */
                    random_hypercube(const int num_dim,
                            const uint64_t* ll, const uint64_t* ur)
                    : random_hypercube(gsl_qrng_sobol, num_dim, ll, ur) {
                    }

                    /**
                     * The basic destructor
                     */
                    virtual ~random_hypercube() {
                        if (m_ll) {
                            delete[] m_ll;
                        }
                        if (m_ur) {
                            delete[] m_ur;
                        }
                        if (m_raw_point) {
                            delete[] m_raw_point;
                        }
                        if (m_is_qr_seq_own && m_qr_seq) {
                            gsl_qrng_free(m_qr_seq);
                        }
                        if (m_scale) {
                            delete[] m_scale;
                        }
                        if (m_shift) {
                            delete[] m_shift;
                        }
                    }

                    /**
                     * Allows to get the vector of lower left hypercube borders
                     * @return the vector of lower left hypercube borders
                     */
                    const uint64_t* get_ll() const {
                        return m_ll;
                    }

                    /**
                     * Allows to get the vector of upper right hypercube borders
                     * @return the vector of upper right hypercube borders
                     */
                    const uint64_t* get_ur() const {
                        return m_ur;
                    }

                    /**
                     * Allows to get the number of hypercube dimensions
                     * @return the number of hypercube dimensions
                     */
                    inline int get_num_dim() const {
                        return m_num_dim;
                    }

                    /**
                     * Allow to generate a new random point, not thread safe!
                     * @param dof_ids_abs the array to be filled in with the sample dof ids
                     * @param dof_ids_dbl the array to be filled in with the sample dof ids
                     */
                    void get_random(uint64_t * dof_ids_abs, double * dof_ids_dbl) {
                        //Get a new raw random point in a thread-safe way
                        {
                            lock_guard<mutex> guard(m_mutex);
                            gsl_qrng_get(m_qr_seq, m_raw_point);
                        }

                        //The point is in the (0,1) range, so scale, shift and round
                        for (int dof_idx = 0; dof_idx < m_num_dim; ++dof_idx) {
                            dof_ids_abs[dof_idx] = round(m_raw_point[dof_idx] * m_scale[dof_idx] + m_shift[dof_idx]);
                            dof_ids_dbl[dof_idx] = dof_ids_abs[dof_idx];
                            //LOG("Generated " << dof_idx << " value: "
                            //        << dof_ids_abs[dof_idx] << " from [" << m_ll[dof_idx]
                            //        << ", " << m_ur[dof_idx] << "]");
                        }
                    }

                    /**
                     * Allows to split the hypercube into two pieces.
                     * @param dof_idx the dimension to bisect in
                     * @param lrb the right border of the left interval to be set
                     * @param rlb the left border of the right interval to be set
                     * @param left_cube the new left part of the hypercube
                     * @param right_cube the new right part of the hypercube
                     * @return true if the bisection was successful, otherwise false
                     */
                    void split(const int dof_idx,
                            const uint64_t lrb, const uint64_t rlb, 
                            random_hypercube*&left_cube,
                            random_hypercube*&right_cube) {
                        //Allocate and initialize the left right interval border arrays
                        const size_t data_size = m_num_dim * sizeof (uint64_t);

                        //Make data for the left interval
                        uint64_t* l_ll = new uint64_t[m_num_dim];
                        memcpy(l_ll, m_ll, data_size);
                        uint64_t* l_ur = new uint64_t[m_num_dim];
                        memcpy(l_ur, m_ur, data_size);
                        l_ur[dof_idx] = lrb;

                        //Make data for the right interval
                        uint64_t* r_ll = new uint64_t[m_num_dim];
                        memcpy(r_ll, m_ll, data_size);
                        uint64_t* r_ur = new uint64_t[m_num_dim];
                        memcpy(r_ur, m_ur, data_size);
                        r_ll[dof_idx] = rlb;

                        //Allocate the new random hypercubes
                        left_cube = new random_hypercube(
                                m_qr_seq, false, m_num_dim, l_ll, l_ur);
                        right_cube = new random_hypercube(
                                m_qr_seq, false, m_num_dim, r_ll, r_ur);
                    }
                };
            }
        }
    }
}

#endif /* RANDOM_HYPERCUBE_HH */

