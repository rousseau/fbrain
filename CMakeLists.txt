PROJECT(fbrain)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#---- Options --------------------------------------------------------------------------

# Parallel computing
OPTION(USE_OMP "Use OpenMP" ON)
IF(USE_OMP)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp  ")
ENDIF(USE_OMP)

# Testing
OPTION(BUILD_TESTING "Build testing applications" OFF)

# Doxygen documentation
INCLUDE(${CMAKE_ROOT}/Modules/Documentation.cmake)
OPTION(BUILD_DOCUMENTATION "Build the doxygen documentation" OFF)

#---- Dependancies ---------------------------------------------------------------------

# Find ITK
MESSAGE(STATUS "Check for ITK package")
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

# Find VTK package
MESSAGE(STATUS "Check for VTK package")
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# Find ANN
FIND_PATH(ANN_INCLUDE_DIR ANN.h /usr/include/ANN)
IF(ANN_INCLUDE_DIR-NOTFOUND)
  MESSAGE(FATAL_ERROR "ANN header not found !")
ENDIF(ANN_INCLUDE_DIR-NOTFOUND)

SET(ANN_NAMES ${ANN_NAMES} ann)
FIND_LIBRARY(ANN_LIBRARY NAMES ${ANN_NAMES} )
IF(ANN_LIBRARY-NOTFOUND)
  MESSAGE(FATAL_ERROR "ANN library not found !")
ENDIF(ANN_LIBRARY-NOTFOUND)

# Find TCLAP
MESSAGE(STATUS "Check for TCLAP")
FIND_PATH(TCLAP_DIRECTORY CmdLine.h /usr/include/tclap)
IF(TCLAP_DIRECTORY-NOTFOUND)
  MESSAGE(FATAL_ERROR "TCLAP not found !")
ENDIF(TCLAP_DIRECTORY-NOTFOUND)

#---- Operations -----------------------------------------------------------------------

# Testing
IF(BUILD_TESTING)
    ENABLE_TESTING()
    SUBDIRS(
        Tests
    )
ENDIF(BUILD_TESTING)

# Documentation
IF(BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
    IF(DOXYGEN_FOUND)
        CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Documentation/Doxygen/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Documentation/Doxygen/Doxyfile @ONLY)
        ADD_CUSTOM_TARGET(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Documentation/Doxygen/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT
            "Generating API documentation with Doxygen" VERBATIM)
    ENDIF(DOXYGEN_FOUND)
ENDIF(BUILD_DOCUMENTATION)

# Executables
SUBDIRS(
    Code
    Applications
    Utilities
)
