#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release_Linux
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/scots2dll.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -m64 -O3
CXXFLAGS=-m64 -m64 -O3

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/usr/local/lib -lcudd /usr/local/lib/libgsl.a

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/libSCOTS2DLL.${CND_DLIB_EXT}

${CND_DISTDIR}/libSCOTS2DLL.${CND_DLIB_EXT}: /usr/local/lib/libgsl.a

${CND_DISTDIR}/libSCOTS2DLL.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}
	${LINK.cc} -o ${CND_DISTDIR}/libSCOTS2DLL.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -m64 -O3 -shared -fPIC

${OBJECTDIR}/src/scots2dll.o: src/scots2dll.cc nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Werror -I/usr/local/jdk1.8.0_171/include/linux/ -I/usr/local/jdk1.8.0_171/include/ -I/usr/local/include -Iinc -I../SCOTS2JNI/target/jni -I../SCOTS2C/ext/SCOTSv2.0/src -I../SCOTS2C/ext/SCOTSv2.0/utils -I../SCOTS2C/src/optdet -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/scots2dll.o src/scots2dll.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
