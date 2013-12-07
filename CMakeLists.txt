#
#GoGoRPIM - A RPIM implementation and 2-D Eletromagnetic simulator
#    Copyright (C) 2012  Péricles Lopes Machado (LANE-UFPA)

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#
# Lane Maxwell Meshless Simulator
# Copyright (c) 2012 Péricles Lopes Machado
#
# CMakeLists - CMake recipe
#

cmake_minimum_required(VERSION 2.6)

# Project details
# ---------------------------------------------------------------------------
project(Lane_Maxwell)

# Include directories
# ---------------------------------------------------------------------------
include_directories(.)

# Packages
# ---------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} .)
find_package(PkgConfig)

# Some variables
# ---------------------------------------------------------------------------
set(LANE_MAXWELL_CFLAGS)
set(LANE_MAXWELL_LDFLAGS)
set(LANE_MAXWELL_LIBRARIES)
set(LANE_MAXWELL_INCLUDE_DIRS)

# 3rd-party libraries
# ---------------------------------------------------------------------------
#include(CMakeLists.depends)

# Environment checks
# ---------------------------------------------------------------------------
if (WIN32)
	add_definitions(-DLANE_MAXWELL_WIN32)
endif (WIN32)

if (APPLE)
	add_definitions(-DLANE_MAXWELL_APPLE)
endif (APPLE)

if (CYGWIN)
	add_definitions(-DLANE_MAXWELL_CYGWIN)
endif (CYGWIN)

if (MSVC)
	add_definitions(-DLANE_MAXWELL_MSVC)
	add_definitions(-DLANE_MAXWELL_MSVC_VERSION=${MSVC_VERSION})
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)
endif (MSVC)

# Build environment
# ---------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Debug)
	message(STATUS "No build specified, defaulting to Debug. Use -DCMAKE_BUILD_TYPE=Release to change")
endif (NOT CMAKE_BUILD_TYPE)


if (MSVC)
	set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} /W1 /wd4503 /GR- /EHsc-")
	set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} /Zi /D_DEBUG /DCLEVER_DEBUG /MDd")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG /MD")
else (MSVC)
	set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall -ansi -fno-rtti -fno-exceptions")
	set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -fno-inline -ggdb -D_DEBUG -DCLEVER_DEBUG")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")
endif (MSVC)




# Matrix library
# ---------------------------------------------------------------------------

# Static library

set(MATRIX_SRC
	math/matrix.cc
	math/point.cc
	domains/domains.cc
	domains/supportdomain.cc
)

add_library(matrix 
	STATIC 
	${MATRIX_SRC}
)

set_target_properties (matrix
	PROPERTIES
	VERSION 1.0.0
	SOVERSION 1)
	

# Shared library
if (NOT WIN32)
	add_library(matrix_dl 
		SHARED
		${MATRIX_SRC}
	)

	set_target_properties (matrix_dl
		PROPERTIES
		VERSION 1.0.0
		SOVERSION 1)

	install(TARGETS matrix matrix_dl
	  RUNTIME DESTINATION bin
	  LIBRARY DESTINATION lib
	  ARCHIVE DESTINATION lib
	)
endif (NOT WIN32)


# RPIM library
# ---------------------------------------------------------------------------

# Static library

set(RPIM_SRC
	rpim/rpim.cc
)

add_library(rpim
	STATIC 
	${RPIM_SRC}
)

set_target_properties (rpim
	PROPERTIES
	VERSION 1.0.0
	SOVERSION 1)
	

# Shared library
if (NOT WIN32)
	add_library(rpim_dl 
		SHARED
		${RPIM_SRC}
	)

	set_target_properties (rpim_dl
		PROPERTIES
		VERSION 1.0.0
		SOVERSION 1)

	install(TARGETS rpim rpim_dl
	  RUNTIME DESTINATION bin
	  LIBRARY DESTINATION lib
	  ARCHIVE DESTINATION lib
	)
endif (NOT WIN32)

# MAXWELL library
# ---------------------------------------------------------------------------

# Static library

set(MAXWELL_SRC
	maxwell/maxwell.cc
	maxwell/maxwellgendom.cc
	maxwell/tmz.cc
	maxwell/tez.cc
	lane/lane.cc
)

add_library(maxwell
	STATIC 
	${MAXWELL_SRC}
)

set_target_properties (maxwell
	PROPERTIES
	VERSION 1.0.0
	SOVERSION 1)
	

# Shared library
if (NOT WIN32)
	add_library(maxwell_dl 
		SHARED
		${MAXWELL_SRC}
	)

	set_target_properties (maxwell_dl
		PROPERTIES
		VERSION 1.0.0
		SOVERSION 1)

	install(TARGETS maxwell maxwell_dl
	  RUNTIME DESTINATION bin
	  LIBRARY DESTINATION lib
	  ARCHIVE DESTINATION lib
	)

endif (NOT WIN32)

# Tests
# ---------------------------------------------------------------------------
set(TEST_SRC
	test_matrix_1
	test_rpim_1
	test_rpim_2
	test_rpim_3
	test_rpim_4
	test_rpim_5
	test_maxwell_1
	test_maxwell_2
	test_maxwell_3
	test_maxwell_4
     test_maxwell_5
     test_maxwell_6
     test_beta
     tradutor1
)

if (WIN32)
	if (NOT PTHREAD_LIB)
		message(STATUS "It's necessary define the muilti-thread library to Windows."
					" Use -DPTHREAD_LIB=pthreadGCE2, for example.")
	endif (NOT PTHREAD_LIB)
endif(WIN32)

if (WIN32)
	foreach (source ${TEST_SRC})
		add_executable(${source}.exe
			tests/${source}.cc
		)

		target_link_libraries(${source}.exe ${LANE_MAXWELL_LIBRARIES} maxwell matrix rpim ${PTHREAD_LIB})
		
		include_directories(${LANE_MAXWELL_INCLUDE_DIRS})
		install(TARGETS ${source}.exe RUNTIME DESTINATION bin)
	endforeach(source)
else (WIN32)
	
	foreach (source ${TEST_SRC})
		add_executable(${source}
			tests/${source}.cc
		)

		target_link_libraries(${source} ${LANE_MAXWELL_LIBRARIES} maxwell matrix rpim pthread)
		include_directories(${LANE_MAXWELL_INCLUDE_DIRS})
		install(TARGETS ${source} RUNTIME DESTINATION bin)
	endforeach(source)
endif(WIN32)

# Executables
# ---------------------------------------------------------------------------

add_executable(genDom
	maxwell/genDom.cc
)

if(WIN32)
    target_link_libraries(genDom ${LANE_MAXWELL_LIBRARIES} maxwell matrix rpim ${PTHREAD_LIB})
else(WIN32)
    target_link_libraries(genDom ${LANE_MAXWELL_LIBRARIES} maxwell matrix rpim pthread)
endif(WIN32)

include_directories(${LANE_MAXWELL_INCLUDE_DIRS})
install(TARGETS genDom RUNTIME DESTINATION bin)

add_executable(maxwellTMZ
	maxwell/maxwellTMZ.cc
)

if(WIN32)
    target_link_libraries(maxwellTMZ ${LANE_MAXWELL_LIBRARIES} maxwell matrix rpim  ${PTHREAD_LIB})
else(WIN32)
    target_link_libraries(maxwellTMZ ${LANE_MAXWELL_LIBRARIES} maxwell matrix rpim pthread)
endif(WIN32)

include_directories(${LANE_MAXWELL_INCLUDE_DIRS})
install(TARGETS maxwellTMZ RUNTIME DESTINATION bin)


add_executable(maxwellTEZ
	maxwell/maxwellTEZ.cc
)


if(WIN32)
#pthreadGCE2
    target_link_libraries(maxwellTEZ ${LANE_MAXWELL_LIBRARIES} maxwell matrix rpim ${PTHREAD_LIB})
else(WIN32)
    target_link_libraries(maxwellTEZ ${LANE_MAXWELL_LIBRARIES} maxwell matrix rpim pthread)
endif(WIN32)

include_directories(${LANE_MAXWELL_INCLUDE_DIRS})
install(TARGETS maxwellTEZ RUNTIME DESTINATION bin)


# Main executable
# ---------------------------------------------------------------------------
#add_executable(clever
#	file.cc
#)
#target_link_libraries(clever ${CLEVER_LIBRARIES})
#include_directories(${CLEVER_INCLUDE_DIRS})
#INSTALL(PROGRAMS bin/maxwell2d/${source} DESTINATION bin)


