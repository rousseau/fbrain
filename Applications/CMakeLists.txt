PROJECT(Applications)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Tools)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Transformations)



# GENERAL PROCESSING
# ------------------
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Denoising)


ADD_EXECUTABLE(btkNLMDenoising btkNLMDenoising.cxx
              ${fbrain_SOURCE_DIR}/Code/Denoising/btkNLMTool.h)

TARGET_LINK_LIBRARIES(btkNLMDenoising ${ITK_LIBRARIES})

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Denoising)

ADD_EXECUTABLE(btkNLMDenoising4DImage btkNLMDenoising4DImage.cxx
               ${fbrain_SOURCE_DIR}/Code/Denoising/btkNLMTool.h)

TARGET_LINK_LIBRARIES(btkNLMDenoising4DImage ${ITK_LIBRARIES})

INSTALL(TARGETS btkNLMDenoising btkNLMDenoising4DImage DESTINATION bin)


# RECONSTRUCTION
# --------------

INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Reconstruction)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/DTI)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Segmentation)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Transformations)

ADD_EXECUTABLE(btkImageReconstruction btkImageReconstruction.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLowToHighImageResolutionMethod.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceRigidRegistration.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceAffineRegistration.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkResampleImageByInjectionFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceTransform.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkImageIntersectionCalculator.h
${fbrain_SOURCE_DIR}/Code/Tools/btkMacro.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRegistration.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkAffineRegistration.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRigidRegistration.h
 )

TARGET_LINK_LIBRARIES(btkImageReconstruction ${ITK_LIBRARIES})


ADD_EXECUTABLE(btkLabelPropagation btkLabelPropagation.cxx
              ${fbrain_SOURCE_DIR}/Code/Segmentation/btkLabelPropagationTool.h)

TARGET_LINK_LIBRARIES(btkLabelPropagation ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSuperResolution btkSuperResolution.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionAffineImageFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionRigidImageFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceTransform.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkMseCostFunction.h
${fbrain_SOURCE_DIR}/Code/Denoising/btkNLMTool.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionRigidImageFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionAffineImageFilter.h
)

TARGET_LINK_LIBRARIES(btkSuperResolution ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSuperResolutionFilterTest btkSuperResolutionFilterTest.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkSuperResolutionFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkSuperResolutionFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkBiasCorrectionFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkBiasCorrectionFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLowToHighImageResolutionMethod.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionReconstructionFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionIBPFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionSRFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionReconstructionFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionIBPFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionSRFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionSliceBySliceAffineFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionSliceBySliceAffineFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionSliceBySliceEulerFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionSliceBySliceEulerFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrection3DAffineFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrection3DAffineFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrection3DEulerFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrection3DEulerFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkSliceRejectionFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkSliceRejectionFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkPSFEstimationFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkPSFEstimationFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceTransform.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkMseCostFunction.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceRigidRegistration.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceAffineRegistration.h
#${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceTransform.h
#${fbrain_SOURCE_DIR}/Code/Reconstruction/btkMseCostFunction.h
${fbrain_SOURCE_DIR}/Code/Tools/btkFileHelper.h
${fbrain_SOURCE_DIR}/Code/Tools/btkFileHelper.cxx
${fbrain_SOURCE_DIR}/Code/Tools/btkImageHelper.h
${fbrain_SOURCE_DIR}/Code/Tools/btkMacro.h


)

TARGET_LINK_LIBRARIES(btkSuperResolutionFilterTest ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkGroupwiseS2SDistortionCorrection btkGroupwiseS2SDistortionCorrection.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolation.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkNumerical.h
${fbrain_SOURCE_DIR}/Code/InputOutput/btkFileNameTools.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkGroupwiseS2SDistortionCorrection.h
 )

TARGET_LINK_LIBRARIES(btkGroupwiseS2SDistortionCorrection ${ITK_LIBRARIES} ${ANN_LIBRARY})



TARGET_LINK_LIBRARIES(btkGroupwiseS2SDistortionCorrection ${ITK_LIBRARIES} ${ANN_LIBRARY})

ADD_EXECUTABLE(btkRBFInterpolationS2S btkRBFInterpolationS2S.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolation.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.cxx
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkNumerical.h
${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolateImageFunctionS2S.h
${fbrain_SOURCE_DIR}/Code/DTI/btkDiffusionGradientTable.h
${fbrain_SOURCE_DIR}/Code/InputOutput/btkFileNameTools.h
)

TARGET_LINK_LIBRARIES(btkRBFInterpolationS2S ${ITK_LIBRARIES} ${ANN_LIBRARY})



INSTALL(TARGETS btkImageReconstruction
                btkLabelPropagation
                btkSuperResolution
                btkSuperResolutionFilterTest
                btkGroupwiseS2SDistortionCorrection
                btkRBFInterpolationS2S
        DESTINATION bin)


# TRACTOGRAPHY
# ------------

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Tractography)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/InputOutput)

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Maths)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Diffusion)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Tools)

ADD_EXECUTABLE(btkTractography btkTractography.cxx)
TARGET_LINK_LIBRARIES(btkTractography btkToolsLibrary btkMathsLibrary btkDiffusionLibrary btkTractographyLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkODFParticleFilterTractography
        btkODFParticleFilterTractography.cxx
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkTypes.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkCartesianCoordinates.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkCartesianCoordinates.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalCoordinates.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalCoordinates.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkVector.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkVector.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkDirection.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkDirection.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkPoint.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkPoint.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSignalExtractor.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkSignalExtractor.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSignal.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkSignal.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalHarmonics.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkSphericalHarmonics.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSHModel.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkSHModel.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkSHModelEstimator.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkSHModelEstimator.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkAPrioriDensity.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkAPrioriDensity.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkLikelihoodDensity.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkLikelihoodDensity.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkImportanceDensity.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkImportanceDensity.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkParticle.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkParticle.h
        ${fbrain_SOURCE_DIR}/Code/Tractography/btkParticleFilter.cxx ${fbrain_SOURCE_DIR}/Code/Tractography/btkParticleFilter.h
)
TARGET_LINK_LIBRARIES(btkODFParticleFilterTractography ${ITK_LIBRARIES} vtkHybrid)

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/DtiParticleFilteringTractography)

ADD_EXECUTABLE(btkTensorParticleFilterTractography
        btkTensorParticleFilterTractography.cxx
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

TARGET_LINK_LIBRARIES(btkTensorParticleFilterTractography ${ITK_LIBRARIES} vtkHybrid)

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

INSTALL(TARGETS btkODFParticleFilterTractography
                btkODFStreamlineTractography
                btkTensorParticleFilterTractography
                btkTensorStreamlineTractography
        DESTINATION bin)

