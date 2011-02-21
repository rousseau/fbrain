PROJECT(fbrain)

cmake_minimum_required(VERSION 2.6)

#------------------------------------------------------------------------------

OPTION(USE_OMP "Use OpenMP" ON)
IF(USE_OMP)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp  ")
ENDIF(USE_OMP)

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
MESSAGE(STATUS "Check for ANN library")
FIND_LIBRARY(ANN_LIBRARY ANN)
IF(ANN_LIBRARY-NOTFOUND)
    MESSAGE(FATAL_ERROR "ANN library not found !")
ENDIF(ANN_LIBRARY-NOTFOUND)
MESSAGE(STATUS "Check for ANN headers")
FIND_PATH(ANN_INCLUDE_DIR ANN.h /usr/include/ANN)
IF(ANN_LIBRARY-NOTFOUND)
    MESSAGE(FATAL_ERROR "ANN headers not found !")
ENDIF(ANN_LIBRARY-NOTFOUND)

# Find TCLAP
MESSAGE(STATUS "Check for TCLAP")
FIND_PATH(TCLAP_DIRECTORY CmdLine.h /usr/include/tclap)
IF(TCLAP_DIRECTORY-NOTFOUND)
  MESSAGE(FATAL_ERROR "TCLAP not found !")
ENDIF(TCLAP_DIRECTORY-NOTFOUND)

#------------------------------------------------------------------------------

SUBDIRS(
  Utilities
  Applications
)


