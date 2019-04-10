# CMake file for Jetscape analysis.
# Based heavily on JETSCAPE CMake file.

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(libJetscape-Analysis)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "./lib")

#######################################

FILE(GLOB JetscapeAnalysis *.cxx)
set (SOURCES  ${SOURCES}  ${JetscapeAnalysis} )

add_library(JetscapeAnalysis SHARED ${SOURCES})
target_link_libraries(JetscapeAnalysis ${ROOT_LIBRARIES} -L${HEPMC_LIB_DIR} -lHepMC)
