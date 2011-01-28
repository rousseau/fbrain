PROJECT(Applications)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNLMDenoising btkNLMDenoising.cxx )
TARGET_LINK_LIBRARIES(btkNLMDenoising ITKIO)

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Reconstruction)

ADD_EXECUTABLE(btkImageReconstruction btkImageReconstruction.cxx )
TARGET_LINK_LIBRARIES(btkImageReconstruction ITKIO ITKNumerics)

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Tractography)

ADD_EXECUTABLE(btkTractography
        btkTractography.cxx
        btkTypes.h
        btkCartesianCoordinates.cxx btkCartesianCoordinates.h
        btkSphericalCoordinates.cxx btkSphericalCoordinates.h
        btkVector.cxx btkVector.h
        btkDirection.cxx btkDirection.h
        btkPoint.cxx btkPoint.h
        btkSignalExtractor.cxx btkSignalExtractor.h
        btkSignal.cxx btkSignal.h
        btkSphericalHarmonics.cxx btkSphericalHarmonics.h
        btkSHModel.cxx btkSHModel.h
        btkSHModelEstimator.cxx btkSHModelEstimator.h
        btkDensity.h
        btkSimDensity.h
        btkSHModelDensity.cxx btkSHModelDensity.h
        btkVonMisesFisherDensity.cxx btkVonMisesFisherDensity.h
        btkNormalDensity.cxx btkNormalDensity.h
        btkInitialDensity.cxx btkInitialDensity.h
        btkAPrioriDensity.cxx btkAPrioriDensity.h
        btkLikelihoodDensity.cxx btkLikelihoodDensity.h
        btkImportanceDensity.cxx btkImportanceDensity.h
        btkParticle.cxx btkParticle.h
        btkParticleFilter.cxx btkParticleFilter.h
)

TARGET_LINK_LIBRARIES(btkTractography ITKCommon ITKIO vtkHybrid)

#------------------------------------------------------------------------------

INSTALL(TARGETS btkNLMDenoising
				btkImageReconstruction
                btkTractography
        DESTINATION bin)

