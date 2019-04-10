# CMake file for Jetscape analysis.
# Based heavily on JETSCAPE CMake file.
# You may need to set the CMAKE_MODULE_PATH to point to JETSCAPE/cmakemodules

###########################
### Initial Cmake Setup ###
###########################

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(Jetscape-analysis CXX)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
message("source directory:" ${CMAKE_SOURCE_DIR})
message("binary directory:" ${CMAKE_BINARY_DIR})
message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# for mac compliance
cmake_policy(SET CMP0042 NEW)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../JETSCAPE/cmakemodules")

###############################
### Compiler & Linker Flags ###
###############################
message("Compiler and Linker flags ...")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -fPIC -pipe -Wall -std=c++11")
## can turn off some warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder -Wno-unused-variable ")
## can turn on debugging information
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

### system dependence.
### note that APPLE also defines UNIX, hence the elseif to differentiate
if (APPLE)
    message( STATUS "Apple : " ${CMAKE_HOST_SYSTEM})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") 
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-but-set-variable -Wno-parentheses -fext-numeric-literals")
    endif()

    ## can turn off some warnings
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-private-field")
elseif(UNIX) 
    message( STATUS "Linux : " ${CMAKE_HOST_SYSTEM})
    ## can turn off some warnings
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-but-set-variable -Wno-parentheses -fext-numeric-literals")
    endif()

    ## Additional options
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()   

message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
#message(STATUS "LD_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS})

# Tell external library path
LINK_DIRECTORIES(/usr/local/lib)


#############################################
### Installing Header and Library Files   ###
#############################################
# Install header files
# default install prefix: build directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}"
    CACHE PATH "Install path prefix, prepended onto install directories."
    FORCE)
endif()
message(STATUS "Now: install prefix is ${CMAKE_INSTALL_PREFIX}")

install(
  # source directory
  DIRECTORY
  # our source
  "${CMAKE_SOURCE_DIR}/src/"
  DESTINATION "include" # target directory
  FILES_MATCHING # install only matched files
  PATTERN "*.h*" # select header files
  )

# Install lib files
install(
  # our libraries
  DIRECTORY
  "${PROJECT_BINARY_DIR}/src/lib/"
  "${PROJECT_BINARY_DIR}/lib/"
  DESTINATION "lib" # target directory
  FILES_MATCHING # install only matched files
  PATTERN "lib*.*" # selects .so, .a, .dylib, ...
  PATTERN "CMakeFiles" EXCLUDE
  )

###########################
### Optional packages   ###
###########################
unset(USE_HEPMC)
message("Looking for HepMC ...")
find_package(HEPMC)
if (${HEPMC_FOUND})
include_directories(${HEPMC_INCLUDE_DIR})
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DUSE_HEPMC" )
Message(STATUS "HepMC Include dir : " ${HEPMC_INCLUDE_DIR})
endif()

message("Looking for ROOT ...")
find_package(ROOT)
if (${ROOT_FOUND})
include_directories(${ROOT_INCLUDES})
Message(STATUS "ROOT Include dir : " ${ROOT_INCLUDES})
Message(STATUS "ROOT Libraries : " ${ROOT_LIBRARIES})
endif()

###########################
### Framework Includes  ###
###########################
message ("Include Directories ...")
include_directories(./src/ )

### include_directories will be expanded as packages are found

###################################################
### Some additional settings for subdirectories ###
###################################################
add_subdirectory(./src)

###########################
### Binary location     ###
###########################
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/..)

###########################
### Executables         ###
###########################

### Run Jetscape
add_executable(runJetscapeAnalysis ./src/runJetscapeAnalysis.C)
target_link_libraries(runJetscapeAnalysis JetscapeAnalysis ${ROOT_LIBRARIES} -L${HEPMC_LIB_DIR} -lHepMC)
