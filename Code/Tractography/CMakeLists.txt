#
#    Name: ./CMakeLists.txt
# Project: BTK Tractography
#


CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
PROJECT(btkTractography)

# -----------------------------------------------------------------------------

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#------------------------------------------------------------------------------

# Find ITK package
MESSAGE(STATUS "Check for ITK package")
FIND_PACKAGE(ITK REQUIRED)
MESSAGE(STATUS "ITK package found")
INCLUDE(${ITK_USE_FILE})

# Find VTK package
MESSAGE(STATUS "Check for VTK package")
FIND_PACKAGE(VTK REQUIRED)
MESSAGE(STATUS "VTK package found")
INCLUDE(${VTK_USE_FILE})

# Find TCLAP
MESSAGE(STATUS "Check for TCLAP")
FIND_PATH(TCLAP_DIRECTORY CmdLine.h)

IF(TCLAP_DIRECTORY-NOTFOUND)
	MESSAGE(FATAL_ERROR "TCLAP not found !")
ENDIF(TCLAP_DIRECTORY-NOTFOUND)

#------------------------------------------------------------------------------

# output path
SET(EXECUTABLE_OUTPUT_PATH
    ${btkTractography_BINARY_DIR}/bin
    CACHE INTERNAL
    "Single output directory for building all executables.")
SET(LIBRARY_OUTPUT_PATH
    ${btkTractography_BINARY_DIR}/lib
    CACHE INTERNAL
    "Single output directory for building all libraries.")

#------------------------------------------------------------------------------

# Adding sources directories
#INCLUDE_DIRECTORIES(tclap)
INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})
INCLUDE_DIRECTORIES(btkCommon)
INCLUDE_DIRECTORIES(btkSignal)
INCLUDE_DIRECTORIES(btkModel)
INCLUDE_DIRECTORIES(btkTractography)


# ---- Libraries --------------------------------------------------------------

# btkCommon library
ADD_LIBRARY(btkCommon STATIC
        btkCommon/btkTypes.h
        btkCommon/btkCartesianCoordinates.cxx btkCommon/btkCartesianCoordinates.h
        btkCommon/btkSphericalCoordinates.cxx btkCommon/btkSphericalCoordinates.h
        btkCommon/btkVector.cxx btkCommon/btkVector.h
        btkCommon/btkDirection.cxx btkCommon/btkDirection.h
        btkCommon/btkPoint.cxx btkCommon/btkPoint.h
)


# btkSignal library
ADD_LIBRARY(btkSignal STATIC
        btkSignal/btkSignalExtractor.cxx btkSignal/btkSignalExtractor.h
        btkSignal/btkSignal.cxx btkSignal/btkSignal.h)

TARGET_LINK_LIBRARIES(btkSignal btkCommon)
TARGET_LINK_LIBRARIES(btkSignal ITKCommon ITKIO)


# btkModel library
ADD_LIBRARY(btkModel STATIC
        btkModel/btkSphericalHarmonics.cxx btkModel/btkSphericalHarmonics.h
        btkModel/btkSHModel.cxx btkModel/btkSHModel.h
        btkModel/btkSHModelEstimator.cxx btkModel/btkSHModelEstimator.h)

TARGET_LINK_LIBRARIES(btkModel btkCommon)
TARGET_LINK_LIBRARIES(btkModel ITKCommon ITKIO)


# btkParticleFilter library
ADD_LIBRARY(btkParticleFilter STATIC
        btkTractography/btkDensity.h
        btkTractography/btkSimDensity.h
        btkTractography/btkSHModelDensity.cxx btkTractography/btkSHModelDensity.h
        btkTractography/btkVonMisesFisherDensity.cxx btkTractography/btkVonMisesFisherDensity.h
        btkTractography/btkNormalDensity.cxx btkTractography/btkNormalDensity.h
        btkTractography/btkInitialDensity.cxx btkTractography/btkInitialDensity.h
        btkTractography/btkAPrioriDensity.cxx btkTractography/btkAPrioriDensity.h
        btkTractography/btkLikelihoodDensity.cxx btkTractography/btkLikelihoodDensity.h
        btkTractography/btkImportanceDensity.cxx btkTractography/btkImportanceDensity.h
        btkTractography/btkParticle.cxx btkTractography/btkParticle.h
        btkTractography/btkParticleFilter.cxx btkTractography/btkParticleFilter.h)

TARGET_LINK_LIBRARIES(btkParticleFilter vtkHybrid)


# ---- Executables ------------------------------------------------------------

# Signal extraction executable
ADD_EXECUTABLE(btkSignalExtraction
        btkSignalExtraction.cxx)

TARGET_LINK_LIBRARIES(btkSignalExtraction btkCommon btkSignal)


# Model estimation executable
ADD_EXECUTABLE(btkModelEstimation
        btkModelEstimation.cxx)

TARGET_LINK_LIBRARIES(btkModelEstimation btkCommon btkModel)


# Tractography executable
ADD_EXECUTABLE(btkParticleFiltering
        btkParticleFiltering.cxx)

TARGET_LINK_LIBRARIES(btkParticleFiltering btkCommon btkSignal btkModel btkParticleFilter)


# Tractography executable
ADD_EXECUTABLE(btkTractography
        btkTractography.cxx
)

TARGET_LINK_LIBRARIES(btkTractography btkCommon btkSignal btkModel btkParticleFilter)
TARGET_LINK_LIBRARIES(btkTractography ITKCommon ITKIO)
TARGET_LINK_LIBRARIES(btkTractography vtkHybrid)


# ---- Installation directories -----------------------------------------------

INSTALL(TARGETS btkSignalExtraction
                btkModelEstimation
                btkParticleFiltering
                btkTractography
        DESTINATION bin)
