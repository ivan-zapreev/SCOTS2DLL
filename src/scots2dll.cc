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

#include "ctrl_wrapper.hh"
#include "ftn_computer.hh"

using namespace std;
using namespace tud::ctrl::scots::jni;

//Stores the fitness class name
static const char * const FITNESS_CLASS_NAME = "nl/tudelft/dcsc/scots2sr/sr/ScaledFitness";
//Stores the fitness computer
static ftn_computer * m_p_ftn_comp = NULL;

JNIEXPORT jint JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_load
(JNIEnv * env, jclass, jstring file_name) {
    /*Re-set the data*/
    if (m_p_ftn_comp) {
        delete m_p_ftn_comp;
    }
    //Instantiate the new fitness computer
    m_p_ftn_comp = new ftn_computer();

    //Load the controller
    const char *file_name_c = env->GetStringUTFChars(file_name, 0);
    const int result = m_p_ftn_comp->load(env, file_name_c);
    env->ReleaseStringUTFChars(file_name, file_name_c);

    return result;
}

JNIEXPORT void JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_configure
(JNIEnv * env, jclass, jobject cfg) {
    (*m_p_ftn_comp) << "Start configuring the object\n" << std::flush;
    m_p_ftn_comp->configure(env, cfg);
}

static jobject compute_fitness(JNIEnv * env, jstring pclass_name, jint ipt_dof_idx) {
    double exact_ftn = 0.0, req_ftn = 0.0, scale = 1.0, shift = 0.0;
    const char *pname_c = env->GetStringUTFChars(pclass_name, 0);
    string pn_local(pname_c);
    jclass ind_class = env->FindClass(pname_c);
    env->ReleaseStringUTFChars(pclass_name, pname_c);
    if (ind_class != NULL) {
        ctrl_wrapper wrapper(env, ind_class, pn_local, ipt_dof_idx);
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
    return env->NewObject(m_ftn_cls, m_ftn_con_id, exact_ftn, req_ftn, scale, shift);
}

JNIEXPORT jobject JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_compute_1fitness
(JNIEnv * env, jclass, jstring person_name, jint ipt_dof_idx) {
    return compute_fitness(env, person_name, ipt_dof_idx);
}

JNIEXPORT void JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_start_1unfit_1export
(JNIEnv * env, jclass) {
    m_p_ftn_comp->start_unfit_export(env);
}

JNIEXPORT void JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_export_1unfit_1points
(JNIEnv * env, jclass, jstring pclass_name, jint ipt_dof_idx) {
    const char *pname_c = env->GetStringUTFChars(pclass_name, 0);
    string pn_local(pname_c);
    jclass ind_class = env->FindClass(pname_c);
    env->ReleaseStringUTFChars(pclass_name, pname_c);
    if (ind_class != NULL) {
        ctrl_wrapper wrapper(env, ind_class, pn_local, ipt_dof_idx);
        m_p_ftn_comp->compute_unfit_points(env, wrapper);
    } else {
        (void) throwException(env, ClassNotFoundException,
                "The requested individual is not found!");
    }
}

JNIEXPORT void JNICALL Java_nl_tudelft_dcsc_scots2jni_Scots2JNI_finish_1unfit_1export
(JNIEnv * env, jclass, jstring file_name) {
    const char *file_name_c = env->GetStringUTFChars(file_name, 0);
    m_p_ftn_comp->finish_unfit_export(env, file_name_c);
    env->ReleaseStringUTFChars(file_name, file_name_c);
}

