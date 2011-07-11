PROJECT(Applications)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNLMDenoising btkNLMDenoising.cxx )
TARGET_LINK_LIBRARIES(btkNLMDenoising ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNLMDenoising4DImage btkNLMDenoising4DImage.cxx )
TARGET_LINK_LIBRARIES(btkNLMDenoising4DImage ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Reconstruction)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/DTI)

ADD_EXECUTABLE(btkImageReconstruction btkImageReconstruction.cxx )
TARGET_LINK_LIBRARIES(btkImageReconstruction ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSuperResolution btkSuperResolution.cxx )
TARGET_LINK_LIBRARIES(btkSuperResolution ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkGroupwiseS2SDistortionCorrection btkGroupwiseS2SDistortionCorrection.cxx )
TARGET_LINK_LIBRARIES(btkGroupwiseS2SDistortionCorrection ${ITK_LIBRARIES} ${ANN_LIBRARY})

ADD_EXECUTABLE(btkRBFInterpolationS2S btkRBFInterpolationS2S.cxx )
TARGET_LINK_LIBRARIES(btkRBFInterpolationS2S ${ITK_LIBRARIES} ${ANN_LIBRARY})

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Tractography)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/InputOutput)

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
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkAPrioriDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkLikelihoodDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkImportanceDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkParticle.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkParticleFilter.cxx
)
TARGET_LINK_LIBRARIES(btkTractography ${ITK_LIBRARIES} vtkHybrid)

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography)

ADD_EXECUTABLE(btkDtiParticleFilteringTractography
        btkDtiParticleFilteringTractography.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkCartesianCoordinates.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalCoordinates.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkVector.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkDirection.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkPoint.cxx
        ${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography/btkDTFPSignalExtractor.cxx
        ${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography/btkDTFPSignal.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalHarmonics.cxx
        ${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography/btkDTFPAPrioriDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography/btkDTFPLikelihoodDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography/btkDTFPImportanceDensity.cxx
        ${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography/btkDTFPParticle.cxx
        ${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography/btkDTFPParticleFilter.cxx
)

TARGET_LINK_LIBRARIES(btkDtiParticleFilteringTractography ${ITK_LIBRARIES} vtkHybrid)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkLabelPropagation btkLabelPropagation.cxx )
TARGET_LINK_LIBRARIES(btkLabelPropagation ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Tractography)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/InputOutput)

ADD_EXECUTABLE(btkODFStreamlineTractography
        btkODFStreamlineTractography.cxx
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
)
TARGET_LINK_LIBRARIES(btkODFStreamlineTractography ${ITK_LIBRARIES} vtkHybrid)

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography)

ADD_EXECUTABLE(btkTensorStreamlineTractography
        btkTensorStreamlineTractography.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkCartesianCoordinates.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalCoordinates.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkVector.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkDirection.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkPoint.cxx
        ${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography/btkDTFPSignalExtractor.cxx
        ${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography/btkDTFPSignal.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalHarmonics.cxx
)

TARGET_LINK_LIBRARIES(btkTensorStreamlineTractography ${ITK_LIBRARIES} vtkHybrid)

#------------------------------------------------------------------------------

INSTALL(TARGETS btkNLMDenoising
                btkNLMDenoising4DImage
                btkImageReconstruction
                btkGroupwiseS2SDistortionCorrection
                btkRBFInterpolationS2S
                btkTractography
                btkLabelPropagation
                btkDtiParticleFilteringTractography
                btkODFStreamlineTractography
                btkTensorStreamlineTractography
        DESTINATION bin)

