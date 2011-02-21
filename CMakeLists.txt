PROJECT(fbrain)

cmake_minimum_required(VERSION 2.6)
OPTION(USE_OMP "Use OpenMP" ON)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp  ")

#------------------------------------------------------------------------------

# Find ITK
MESSAGE(STATUS "Check for ITK package")
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

# Find VTK package
MESSAGE(STATUS "Check for VTK package")
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# Find ANN
MESSAGE(STATUS "Check for ANN package")
FIND_PACKAGE(ANN REQUIRED)
IF(ANN_FOUND)
ELSE(ANN_FOUND)
  MESSAGE(FATAL_ERROR
          "ANN not found. Please perform installation.")
ENDIF(ANN_FOUND)

# Find TCLAP
MESSAGE(STATUS "Check for TCLAP")
FIND_PATH(TCLAP_DIRECTORY CmdLine.h /usr/include/tclap)

IF(TCLAP_DIRECTORY-NOTFOUND)
  MESSAGE(FATAL_ERROR "TCLAP not found !")
ENDIF(TCLAP_DIRECTORY-NOTFOUND)

#------------------------------------------------------------------------------

#SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/bin")

SUBDIRS(
  Utilities
  Applications
)


