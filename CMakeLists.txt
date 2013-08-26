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

# Python script support
OPTION(ENABLE_SCRIPT "Enable the use of python script in BTK." OFF)

# Build Experimental programs
OPTION(BUILD_EXPERIMENTAL "Enable the use of experimental programs in BTK " OFF)

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

# Find CPPunit (Optionnal)
IF(BUILD_TESTING)
    MESSAGE(STATUS "Check for CppUnit")

    FIND_PATH(CPPUNIT_INCLUDE_DIR NAMES TestSuite.h PATHS /usr/include/cppunit /usr/local/include/cppunit ~/.local/include/cppunit)
    IF(CPPUNIT_INCLUDE_DIR-NOTFOUND)
        MESSAGE(FATAL_ERROR "CppUnit include directory not found !")
    ELSE(CPPUNIT_INCLUDE_DIR-NOTFOUND)
        INCLUDE_DIRECTORIES(${CPPUNIT_INCLUDE_DIR})
    ENDIF(CPPUNIT_INCLUDE_DIR-NOTFOUND)

    FIND_PATH(CPPUNIT_LIBRARY NAMES libcppunit.a PATHS /usr/lib /usr/local/lib ~/.local/lib)
    IF(CPPUNIT_LIBRARY-NOTFOUND)
        MESSAGE(FATAL_ERROR "CppUnit library not found !")
    ELSE(CPPUNIT_LIBRARY-NOTFOUND)
        SET(CPPUNIT_LIBRARY ${CPPUNIT_LIBRARY}/libcppunit.a)
    ENDIF(CPPUNIT_LIBRARY-NOTFOUND)
ENDIF(BUILD_TESTING)

#---- Operations -----------------------------------------------------------------------

# Testing
IF(BUILD_TESTING)
    ENABLE_TESTING()
    SUBDIRS(
        Tests
    )
ENDIF(BUILD_TESTING)

# Beta
IF(BUILD_BETA)
    SUBDIRS(
        Experimental
    )
ENDIF(BUILD_BETA)


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

# Uninstall command
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)


# Executables
SUBDIRS(
    Code
    Applications
    Scripts
    Utilities
    Viewer
)
