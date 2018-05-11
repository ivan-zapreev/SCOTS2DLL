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

using namespace std;
using namespace tud::ctrl::scots::jni;

//Stores the extended fitness class name
static const char * const EXT_FITNESS_CLASS_NAME = "nl/tudelft/dcsc/scots2sr/sr/ExtendedFitness";

//Stores the scaled fitness class name
static const char * const SC_FITNESS_CLASS_NAME = "nl/tudelft/dcsc/scots2sr/sr/ScaledFitness";

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

JNIEXPORT void JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_configure
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

static string get_pclass_by_name(JNIEnv * env, jstring pclass_name, jclass & ind_class) {
    const char *pname_c = env->GetStringUTFChars(pclass_name, 0);
    string pn_local(pname_c);
    ind_class = env->FindClass(pname_c);
    env->ReleaseStringUTFChars(pclass_name, pname_c);
    if (ind_class == NULL) {
        (void) throwException(env, ClassNotFoundException,
                "The requested individual is not found!");
    }
    return pn_local;
}

/**
 * Allows to compute plain fitness with scaling
 * @param env the JNI environment
 * @param wrap the individual clas wrapper
 * @return the fitness object
 */
static jobject compute_fitness_sc(
        JNIEnv * env, ctrl_wrapper & wrap) {
    //Create java arrays for storing scales and shifts
    const int is_dim = wrap.get_is_dim();
    jdoubleArray jscale = env->NewDoubleArray(is_dim);
    jdoubleArray jshift = env->NewDoubleArray(is_dim);

    //Access the actual array elements to fill in
    double * scale = env->GetDoubleArrayElements(jscale, 0);
    double * shift = env->GetDoubleArrayElements(jshift, 0);

    //Compute the scaled fitness, with the scales and shifts
    double act_ftn = 0.0, ext_ftn = 0.0;
    m_p_ftn_comp->compute(env, wrap, act_ftn, ext_ftn, scale, shift);

    //Release the elements to be copied to java
    env->ReleaseDoubleArrayElements(jscale, scale, 0);
    env->ReleaseDoubleArrayElements(jshift, shift, 0);

    //Find the scaled fitness container class
    jclass scf_cls = env->FindClass(SC_FITNESS_CLASS_NAME);
    //Find the scaled fitness container class constructor
    jmethodID scf_con = env->GetMethodID(scf_cls, "<init>", "(DD[D[D)V");
    //Instantiate the result
    return env->NewObject(scf_cls, scf_con, act_ftn, ext_ftn, jscale, jshift);
}

/**
 * Allows to compute plain fitness without scaling
 * @param env the JNI environment
 * @param wrap the individual clas wrapper
 * @return the fitness object
 */
static jobject compute_fitness_pl(
        JNIEnv * env, ctrl_wrapper & wrap) {
    //Compute the scaled fitness
    double act_ftn = 0.0, ext_ftn = 0.0;
    m_p_ftn_comp->compute(env, wrap, act_ftn, ext_ftn);

    //Find the extended fitness container class
    jclass exf_cls = env->FindClass(EXT_FITNESS_CLASS_NAME);
    //Find the extended fitness container class constructor
    jmethodID exf_con = env->GetMethodID(exf_cls, "<init>", "(DD)V");
    //Instantiate the result
    return env->NewObject(exf_cls, exf_con, act_ftn, ext_ftn);
}

JNIEXPORT jobject JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_compute_1fitness
(JNIEnv * env, jclass, jstring pclass_name) {
    jobject result = NULL;
    if (m_p_data != NULL) {

        //Obtain the individual's class and name
        jclass ind_class = NULL;
        string pn_local = get_pclass_by_name(env, pclass_name, ind_class);
        const int ss_dim = m_p_data->get_config().m_num_ss_dim;
        const int is_dim = m_p_data->get_config().m_num_is_dim;

        //Instantiate the controller's wrapper
        ctrl_wrapper wrap(env, ind_class, pn_local, ss_dim, is_dim);

        //Compute fitness and the resulting object
        double ex_ftn = 0.0, req_ftn = 0.0;
        if (m_p_data->get_config().m_is_scale) {
            result = compute_fitness_sc(env, wrap);
        } else {
            result = compute_fitness_pl(env, wrap);
        }
    } else {
        (void) throwException(env, IllegalStateException,
                "The controller is not loaded!");
    }

    return result;
}

JNIEXPORT void JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_start_1unfit_1export
(JNIEnv * env, jclass) {
    m_p_unf_exp = new unfit_exporter(*m_p_data);
    m_p_unf_exp->start_unfit_export(env);
}

JNIEXPORT jdouble JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_export_1unfit_1points
(JNIEnv * env, jclass, jstring pclass_name) {
    double result = 0.0;

    if (m_p_data != NULL) {
        //Obtain the individual's class and name
        jclass ind_class = NULL;
        string pn_local = get_pclass_by_name(env, pclass_name, ind_class);

        //Create wrapper class
        ctrl_wrapper wrapper(env, ind_class, pn_local,
                m_p_data->get_config().m_num_ss_dim,
                m_p_data->get_config().m_num_is_dim);

        //Export unfit points
        result = m_p_unf_exp->compute_unfit_points(env, wrapper);
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

