cmake_minimum_required(VERSION 2.8.10)
project(lp-engine)
set( VERSION 0.1 ) 

if( NOT CMAKE_BUILD_TYPE )
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW) # MACOSX_RPATH
endif()


# Include local cmake modules directory.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")
include(install_package)

# Default flags.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall") #-DBOOST_NO_CXX11_NUMERIC_LIMITS")

IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
ENDIF()

option(BUILD_APPLICATIONS "Build Applications" ON)

find_package(Eigen3 REQUIRED)

set( _INCLUDE_DIRS
        ${CMAKE_BINARY_DIR}
        ${EIGEN3_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        )
set( _LIBRARIES
#        ${EIGEN3_LIBRARIES}
        )

include_directories( ${_INCLUDE_DIRS} )

set( HDRS
        include/lp-engine/lpengine.h
        include/lp-engine/lpproblem.h
	include/lp-engine/lpdictionary.h
        )
set( SRCS
        src/lpengine.cpp
	src/lpproblem.cpp
	src/lpdictionary.cpp
        )

if(BUILD_APPLICATIONS)
  add_subdirectory( tests )
endif()

set( LIBRARY_NAME ${PROJECT_NAME} )
string( TOLOWER ${LIBRARY_NAME} LIBRARY_NAME )

add_library( ${LIBRARY_NAME} ${SRCS} ${HDRS} )
target_link_libraries( ${LIBRARY_NAME} ${_LIBRARIES} )

install_package(
  PKG_NAME ${PROJECT_NAME}
  LIB_NAME ${LIBRARY_NAME}
  VERSION ${VERSION}
  DESCRIPTION "LP problem solver"
  INSTALL_INCLUDE_DIR true
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  INCLUDE_DIRS ${_INCLUDE_DIRS}
  LINK_LIBS ${_LIBRARIES}
)

