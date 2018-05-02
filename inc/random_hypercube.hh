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
                     * Allows to bisect the hypercube into two pieces
                     * @param dim_idx the dimension to bisect in
                     * @param left_cube the new left part of the hypercube
                     * @param right_cube the new right part of the hypercube
                     * @return true if the bisection was successful, otherwise false
                     */
                    bool bisect(const int dim_idx,
                            random_hypercube*&left_cube,
                            random_hypercube*&right_cube) {
                        bool result = true;
                        //Get the values
                        const uint64_t left_point = m_ll[dim_idx];
                        const uint64_t right_point = m_ur[dim_idx];
                        const uint64_t num_its = (right_point - left_point);
                        if (num_its >= 3) {
                            //Compute the mid point, is always rounded down
                            const uint64_t mid_point = (right_point + left_point) / 2;
                            //Compute the two new interval bounds
                            uint64_t rb_l = mid_point, lb_r = mid_point;
                            //Check on the number of intervals
                            if (num_its % 2) {
                                //The number of intervals is odd, so the mid point
                                //is  rounded down from the center, thus the left
                                //border of the right interval is to be on the
                                //right of the center
                                lb_r += 1;
                            } else {
                                //The number of intervals is even, so there are
                                //two possibilities to make the bisection, so
                                //make a uniform random choice between them.
                                if (bool_gen()) {
                                    rb_l -= 1;
                                } else {
                                    lb_r += 1;
                                }
                            }

                            //Allocate and initialize the left right interval border arrays
                            const size_t data_size = m_num_dim * sizeof (uint64_t);
                            
                            uint64_t* ll_l = new uint64_t[m_num_dim];
                            memcpy(ll_l, m_ll, data_size);
                            uint64_t* ur_l = new uint64_t[m_num_dim];
                            memcpy(ur_l, m_ur, data_size);
                            ur_l[dim_idx] = rb_l;
                            
                            uint64_t* ll_r = new uint64_t[m_num_dim];
                            memcpy(ll_r, m_ll, data_size);
                            uint64_t* ur_r = new uint64_t[m_num_dim];
                            memcpy(ur_r, m_ur, data_size);
                            ll_r[dim_idx] = lb_r;

                            //Allocate the new random hypercubes
                            left_cube = new random_hypercube(
                                    m_qr_seq, false, m_num_dim, ll_l, ur_l);
                            right_cube = new random_hypercube(
                                    m_qr_seq, false, m_num_dim, ll_r, ur_r);
                        } else {
                            result = false;
                        }
                        return result;
                    }
                };
            }
        }
    }
}

#endif /* RANDOM_HYPERCUBE_HH */

