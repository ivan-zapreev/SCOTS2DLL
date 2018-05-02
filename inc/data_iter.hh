/*
 * Copyright (C) 2018 Dr. Ivan S. Zapreev <ivan.zapreev@gmail.com>
 *
 *  Visit my Linked-in profile:
 *     https://nl.linkedin.com/in/zapreevis
 *  Visit my GitHub:
 *     https://github.com/ivan-zapreev
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * File:   data_iter.hh
 * Author: Dr. Ivan S. Zapreev <ivan.zapreev@gmail.com>
 *
 * Created on April 27, 2018, 3:42 PM
 */

#ifndef DATA_ITER_HH
#define DATA_ITER_HH

namespace tud {
    namespace ctrl {
        namespace scots {
            namespace jni {

                /**
                 * This class allows to iterate through the plain controller data.
                 */
                class data_iter {
                protected:
                    //Stores the configuration wrapper reference
                    const config_wrapper & m_config;
                    //Stores the reference to the data source
                    const data_source & m_data;

                public:

                    /**
                     * The basic constructor
                     * @param data the data source to be used
                     */
                    data_iter(const data_source & data)
                    : m_data(data), m_config(data.get_config()) {
                    }
                };
            }
        }
    }
}

#endif /* DATA_ITER_HH */

