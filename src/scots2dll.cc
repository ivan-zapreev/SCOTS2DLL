/*
 * File:   scots2dll.cc
 * Author: Dr. Ivan S. Zapreev
 *
 * Visit my Linked-in profile:
 *      <https://nl.linkedin.com/in/zapreevis>
 * Visit my GitHub:
 *      <https://github.com/ivan-zapreev>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
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
 * Created on January 08, 2018, 6:49 PM
 */

#include <jni.h>
#include <stdio.h>

#include <cstdint>
#include <vector>
#include <fstream>
#include <cmath>

#include "nl_tudelft_dcsc_scots2jni_Scots2JNI.h"

#include "info_logger.hh"
#include "data_source.hh"
#include "ctrl_wrapper.hh"
#include "fitness_computer.hh"
#include "unfit_exporter.hh"

#define VOID void

using namespace std;
using namespace tud::ctrl::scots::jni;

//Stores the fitness class name
static const char * const FITNESS_CLASS_NAME = "nl/tudelft/dcsc/scots2sr/sr/ScaledFitness";

//Stores the fitness computer
static fitness_computer * m_p_ftn_comp = NULL;

//Stores the pointer to the data source
static data_source * m_p_data = NULL;

//Stores the unfit points exporter
static unfit_exporter * m_p_unf_exp = NULL;

/**
 * Allows to re-set the source file, i.e. re-create 
 * the data source and fitness computer and load the data.
 * @param env the JNI environment
 * @param file_name_c the controller file name
 * @return the number of controller dimensions
 */
static int re_set_source_file(JNIEnv * env, const char *file_name_c) {
    /*Re-set the data*/
    if (m_p_ftn_comp) {
        delete m_p_ftn_comp;
    }
    /*Re-set the data*/
    if (m_p_data) {
        delete m_p_data;
    }

    //Open logging
    open_info_log(file_name_c);

    //Instantiate the new data source
    m_p_data = new data_source();

    //Load the controller
    return m_p_data->load(env, file_name_c);
}

JNIEXPORT jint JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_load
(JNIEnv * env, jclass, jstring file_name) {
    //Get the file name and open logging
    const char *file_name_c = env->GetStringUTFChars(file_name, 0);

    //Re-set the data
    const int result = re_set_source_file(env, file_name_c);

    //Release the string
    env->ReleaseStringUTFChars(file_name, file_name_c);

    return result;
}

JNIEXPORT jint JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_get_1state_1space_1size
(JNIEnv * env, jclass, jint ss_dim) {
    int result = 0;
    if (m_p_data) {
        result = m_p_data->get_state_space_size(env, ss_dim);
    } else {
        (void) throwException(env, IllegalStateException,
                "The controller is not loaded!");
    }

    return result;
}

JNIEXPORT VOID JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_configure
(JNIEnv * env, jclass, jobject cfg) {
    LOG("Start configuring the object");

    if (m_p_data != NULL) {
        //Configure the data source
        m_p_data->configure(env, cfg);

        //Instantiate the new fitness computer, 
        //once the data source is fully configured
        m_p_ftn_comp = new fitness_computer(*m_p_data);
    } else {
        (void) throwException(env, IllegalStateException,
                "The controller is not loaded!");
    }
    LOG("The object configuration is finished");
}

static jobject compute_fitness(JNIEnv * env, jstring pclass_name, jint ipt_dof_idx) {
    jobject result;
    if (m_p_data != NULL) {
        double exact_ftn = 0.0, req_ftn = 0.0, scale = 1.0, shift = 0.0;
        const char *pname_c = env->GetStringUTFChars(pclass_name, 0);
        string pn_local(pname_c);
        jclass ind_class = env->FindClass(pname_c);
        env->ReleaseStringUTFChars(pclass_name, pname_c);
        if (ind_class != NULL) {
            ctrl_wrapper wrapper(env, ind_class, pn_local,
                    ipt_dof_idx, m_p_data->get_config().m_num_ss_dim);
            m_p_ftn_comp->compute(env, wrapper, exact_ftn, req_ftn, scale, shift);
        } else {
            (void) throwException(env, ClassNotFoundException,
                    "The requested individual is not found!");
        }

        //Find the fitness container class
        jclass m_ftn_cls = env->FindClass(FITNESS_CLASS_NAME);
        //Find the fitness container class constructor
        jmethodID m_ftn_con_id = env->GetMethodID(m_ftn_cls, "<init>", "(DDDD)V");
        //Return the fitness object
        result = env->NewObject(m_ftn_cls, m_ftn_con_id, exact_ftn, req_ftn, scale, shift);
    } else {
        (void) throwException(env, IllegalStateException,
                "The controller is not loaded!");
    }
    return result;
}

JNIEXPORT jobject JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_compute_1fitness
(JNIEnv * env, jclass, jstring person_name, jint ipt_dof_idx) {
    return compute_fitness(env, person_name, ipt_dof_idx);
}

JNIEXPORT void JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_start_1unfit_1export
(JNIEnv * env, jclass) {
    m_p_unf_exp = new unfit_exporter(*m_p_data);
    m_p_unf_exp->start_unfit_export(env);
}

JNIEXPORT jdouble JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_export_1unfit_1points
(JNIEnv * env, jclass, jstring pclass_name, jint ipt_dof_idx) {
    double result = 0.0;

    if (m_p_data != NULL) {
        const char *pname_c = env->GetStringUTFChars(pclass_name, 0);
        string pn_local(pname_c);
        jclass ind_class = env->FindClass(pname_c);
        env->ReleaseStringUTFChars(pclass_name, pname_c);
        if (ind_class != NULL) {
            ctrl_wrapper wrapper(env, ind_class, pn_local,
                    ipt_dof_idx, m_p_data->get_config().m_num_ss_dim);
            result = m_p_unf_exp->compute_unfit_points(env, wrapper);
        } else {
            (void) throwException(env, ClassNotFoundException,
                    "The requested individual is not found!");
        }
    } else {
        (void) throwException(env, IllegalStateException,
                "The controller is not loaded!");
    }

    return result;
}

JNIEXPORT void JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_finish_1unfit_1export
(JNIEnv * env, jclass, jstring file_name) {
    const char *file_name_c = env->GetStringUTFChars(file_name, 0);
    m_p_unf_exp->finish_unfit_export(env, file_name_c);
    env->ReleaseStringUTFChars(file_name, file_name_c);
    delete m_p_unf_exp;
}

