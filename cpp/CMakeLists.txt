
# Specify the minimum CMake version required
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# Set the project name
project(LDMX_ANALYSIS)

# Set the current version
set(LDMX_ANALYSIS_VERSION_MAJOR 1)
set(LDMX_ANALYSIS_VERSION_MINOR 0)

# Set the output directories of the binaries, static and shared libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${LDMX_ANALYSIS_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LDMX_ANALYSIS_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${LDMX_ANALYSIS_BINARY_DIR}/lib")

# Set additional compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Set the path to additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# Find ROOT
find_package(ROOT REQUIRED)

# Include the ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Include the LCIO headers
include_directories("$ENV{LCIO}/include")

set(UTILS_HDR_DIR "${LDMX_ANALYSIS_SOURCE_DIR}/include/utils")
include_directories("${UTILS_HDR_DIR}")

set(ANALYSIS_HDR_DIR "${LDMX_ANALYSIS_SOURCE_DIR}/include/analyses")
include_directories("${ANALYSIS_HDR_DIR}")

add_subdirectory("${LDMX_ANALYSIS_SOURCE_DIR}/src/utils")

# Add the analysis subdirectory to the build
add_subdirectory("${LDMX_ANALYSIS_SOURCE_DIR}/src/analyses")
