PROJECT(Applications)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNLMDenoising btkNLMDenoising.cxx )
TARGET_LINK_LIBRARIES(btkNLMDenoising ITKIO)

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Reconstruction)

ADD_EXECUTABLE(btkImageReconstruction btkImageReconstruction.cxx )
TARGET_LINK_LIBRARIES(btkImageReconstruction ITKIO ITKNumerics)

#ADD_EXECUTABLE(btkGroupwiseS2SDistortionCorrection btkGroupwiseS2SDistortionCorrection.cxx )
#TARGET_LINK_LIBRARIES(btkGroupwiseS2SDistortionCorrection ITKIO ITKNumerics ${ANN_LIBRARY})

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Tractography)

ADD_EXECUTABLE(btkTractography
        btkTractography.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkCartesianCoordinates.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalCoordinates.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkVector.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkDirection.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkPoint.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSignalExtractor.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSignal.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalHarmonics.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSHModel.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSHModelEstimator.cxx
#        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSHModelDensity.cxx
#        ${fbrain_SOURCE_DIR}/Code/Tractography/btkInitialDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkAPrioriDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkLikelihoodDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkImportanceDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkParticle.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkParticleFilter.cxx
)

TARGET_LINK_LIBRARIES(btkTractography ITKCommon ITKIO vtkHybrid)

#------------------------------------------------------------------------------

INSTALL(TARGETS btkNLMDenoising
                btkImageReconstruction
                btkTractography
        DESTINATION bin)

