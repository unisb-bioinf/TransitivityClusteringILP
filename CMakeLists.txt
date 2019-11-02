cmake_minimum_required(VERSION 2.6)
project(TransitivityClusteringILP)

####################################################################################################
# Manual entries (edit this for new release)
####################################################################################################

set(TRANSITIVITY_CLUSTERING_ILP_PACKAGE_VERSION_MAJOR "0")
set(TRANSITIVITY_CLUSTERING_ILP_PACKAGE_VERSION_MINOR "1")
set(TRANSITIVITY_CLUSTERING_ILP_PACKAGE_VERSION_PATCH "0")
set(TRANSITIVITY_CLUSTERING_ILP_PACKAGE_VERSION "${TRANSITIVITY_CLUSTERING_ILP_PACKAGE_VERSION_MAJOR}.${TRANSITIVITY_CLUSTERING_ILP_PACKAGE_VERSION_MINOR}.${TRANSITIVITY_CLUSTERING_ILP_PACKAGE_VERSION_PATCH}")

####################################################################################################
# Set the path to the additional cmake files.
####################################################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

include(CompilerSpecific)

####################################################################################################
# Set default build type to Release
####################################################################################################

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

####################################################################################################
# Build Tests
####################################################################################################

# This call needs to be in the top level file before any
# add_subdirectory calls in order to allow the testing
# framework to pick up all unit test.

# We require the gtest src directory to be set
if(NOT GTEST_SRC_DIR)
    message(STATUS "GTEST_SRC_DIR not set, disabling unit tests")
    SET(ENABLE_TESTING FALSE)
else()
    SET(ENABLE_TESTING TRUE)
    enable_testing()
    add_subdirectory(test)
endif()

####################################################################################################
# Find Boost
####################################################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Find Boost
find_package(Boost 1.55.0 COMPONENTS iostreams filesystem system program_options REQUIRED)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost version 1.55.0 or higher is required.")
endif()

SET(Boost_LIBRARIES
	${Boost_IOSTREAMS_LIBRARY}
	${Boost_FILESYSTEM_LIBRARY}
	${Boost_SYSTEM_LIBRARY}
	${Boost_PROGRAM_OPTIONS_LIBRARY}
)

include_directories(${Boost_INCLUDE_DIRS})

SET(TRANSITIVITY_CLUSTERING_ILP_DEP_LIBRARIES
	${Boost_LIBRARIES}
)

####################################################################################################
# Find Eigen3
####################################################################################################

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

####################################################################################################
# Build library and executables
####################################################################################################

include_directories(
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    "${CMAKE_SOURCE_DIR}"
    "${CMAKE_BINARY_DIR}"
    /opt/cplex1262/cplex/include
    /opt/cplex1262/concert/include
)

link_directories(
        /opt/cplex1262/concert/lib/x86-64_linux/static_pic
        /opt/cplex1262/cplex/lib/x86-64_linux/static_pic
)


set(LINK_FLAGS
        "-pthread"
        )

add_library(clustLib SHARED
    macros.h
    Matrix.h                                                                                                                                                                                                                                       
    AbstractMatrix.cpp
    AbstractMatrix.h
    DenseMatrix.cpp                                                                                                                                                                                                                                
    DenseMatrix.h                                                                                                                                                                                                                                  
    DenseMatrixIterator.h                                                                                                                                                                                                                          
    DenseMatrixReader.cpp                                                                                                                                                                                                                          
    DenseMatrixReader.h                                                                                                                                                                                                                            
    DenseMatrixWriter.cpp                                                                                                                                                                                                                          
    DenseMatrixWriter.h                                                                                                                                                                                                                            
    DenseColumnSubset.cpp                                                                                                                                                                                                                          
    DenseColumnSubset.h                                                                                                                                                                                                                            
    DenseRowSubset.cpp                                                                                                                                                                                                                             
    DenseRowSubset.h
    Distance.h                                                                                                                                                                                                                              
    Exception.cpp
    Exception.h
    MatrixIterator.h                                                                                                                                                                                                                              
    MatrixWriter.cpp                                                                                                                                                                                                                               
    MatrixWriter.h 
    matrixTools.cpp
    matrixTools.h
    Statistic.h
    IterativeClusteringILP.cpp
    IterativeClusteringILP.h
    ClusteringILP.cpp
    ClusteringILP.h
)

COMPILE_FLAGS(clustLib)
target_include_directories(clustLib PUBLIC "${CMAKE_BINARY_DIR}/libraries/" clustLib_INCLUDE)
target_link_libraries(clustLib ${TRANSITIVITY_CLUSTERING_ILP_DEP_LIBRARIES}
        ilocplex
        concert
        cplex
        -pthread)

add_executable(transitivityClusteringILP transitivityClusteringILP.cpp)
COMPILE_FLAGS(transitivityClusteringILP)
target_link_libraries(transitivityClusteringILP clustLib ${TRANSITIVITY_CLUSTERING_ILP_DEP_LIBRARIES}
        ilocplex
        concert
        cplex
        -pthread)

set_target_properties(transitivityClusteringILP PROPERTIES
        INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
	/opt/cplex1262/cplex/include
	/opt/cplex1262/concert/include
)
