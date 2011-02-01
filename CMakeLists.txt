PROJECT(fbrain)

cmake_minimum_required(VERSION 2.6)
OPTION(USE_OMP "Use OpenMP" ON)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp  ")

#------------------------------------------------------------------------------

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# Find VTK package
MESSAGE(STATUS "Check for VTK package")
FIND_PACKAGE(VTK REQUIRED)
MESSAGE(STATUS "VTK package found")
INCLUDE(${VTK_USE_FILE})

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


