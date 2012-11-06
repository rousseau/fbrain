PROJECT(Applications)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INCLUDE_DIRECTORIES(
    ${TCLAP_DIRECTORY}
    ${Code_SOURCE_DIR}/Tools
    ${Code_SOURCE_DIR}/Maths
    ${Code_SOURCE_DIR}/Transformations

)


# ---- Denoising ------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Denoising)

ADD_EXECUTABLE(btkNLMDenoising btkNLMDenoising.cxx
    ${fbrain_SOURCE_DIR}/Code/Denoising/btkNLMTool.h
)
TARGET_LINK_LIBRARIES(btkNLMDenoising ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkNLMDenoising4DImage btkNLMDenoising4DImage.cxx
    ${fbrain_SOURCE_DIR}/Code/Denoising/btkNLMTool.h
)
TARGET_LINK_LIBRARIES(btkNLMDenoising4DImage ${ITK_LIBRARIES})

INSTALL(TARGETS
    btkNLMDenoising
    btkNLMDenoising4DImage
DESTINATION bin)

# ---- Reconstruction -------------------------------------------------------------------

INCLUDE_DIRECTORIES(
    ${ANN_INCLUDE_DIR}
    ${Code_SOURCE_DIR}/Reconstruction
    ${Code_SOURCE_DIR}/Reconstruction/SuperResolution
    ${Code_SOURCE_DIR}/Diffusion
    ${Code_SOURCE_DIR}/Segmentation
    ${Code_SOURCE_DIR}/Transformations
    ${Code_SOURCE_DIR}/Registration
)

ADD_EXECUTABLE(btkImageReconstruction btkImageReconstruction.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLowToHighImageResolutionMethod.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceRigidRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceAffineRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkResampleImageByInjectionFilter.h
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkImageIntersectionCalculator.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkAffineRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkRigidRegistration.h
)
TARGET_LINK_LIBRARIES(btkImageReconstruction btkToolsLibrary ${ITK_LIBRARIES})


ADD_EXECUTABLE(btkLabelPropagation btkLabelPropagation.cxx
    ${fbrain_SOURCE_DIR}/Code/Segmentation/btkLabelPropagationTool.h
)
TARGET_LINK_LIBRARIES(btkLabelPropagation ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSuperResolution btkSuperResolution.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionAffineImageFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionRigidImageFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.h
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkMseCostFunction.h
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
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkMseCostFunction.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceRigidRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceAffineRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkMseCostFunction.h
)

TARGET_LINK_LIBRARIES(btkSuperResolutionFilterTest btkToolsLibrary btkMathsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkGroupwiseS2SDistortionCorrection btkGroupwiseS2SDistortionCorrection.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolation.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkNumerical.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkGroupwiseS2SDistortionCorrection.h
)
TARGET_LINK_LIBRARIES(btkGroupwiseS2SDistortionCorrection btkToolsLibrary ${ITK_LIBRARIES} ${ANN_LIBRARY})


ADD_EXECUTABLE(btk3DSuperResolution btk3DSuperResolution.cxx
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
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkMseCostFunction.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceRigidRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceAffineRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkMseCostFunction.h
)

TARGET_LINK_LIBRARIES(btk3DSuperResolution btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkRBFInterpolationS2S btkRBFInterpolationS2S.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolation.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkNumerical.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolateImageFunctionS2S.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionGradientTable.h
)


TARGET_LINK_LIBRARIES(btkRBFInterpolationS2S btkToolsLibrary ${ITK_LIBRARIES} ${ANN_LIBRARY})

ADD_EXECUTABLE(btkIntersectionRegistration btkIntersectionRegistration.cxx
${fbrain_SOURCE_DIR}/Code/Tools/btkImageHelper.h
${fbrain_SOURCE_DIR}/Code/Registration/btkSlicesIntersectionVNLCostFunction.hxx
${fbrain_SOURCE_DIR}/Code/Registration/btkSlicesIntersectionITKCostFunction.hxx)

TARGET_LINK_LIBRARIES(btkIntersectionRegistration btkToolsLibrary btkMathsLibrary ${ITK_LIBRARIES})

INSTALL(TARGETS btkImageReconstruction
                btkLabelPropagation
                btkSuperResolution
                btkSuperResolutionFilterTest
                btkGroupwiseS2SDistortionCorrection
                btkRBFInterpolationS2S
                btkIntersectionRegistration
        DESTINATION bin)



# ---- Tractography ---------------------------------------------------------------------

INCLUDE_DIRECTORIES(
    ${Code_SOURCE_DIR}/Diffusion
    ${Code_SOURCE_DIR}/Tractography
    ${Code_SOURCE_DIR}/DtiParticleFilteringTractography
)

ADD_EXECUTABLE(btkTractography btkTractography.cxx)
TARGET_LINK_LIBRARIES(btkTractography btkToolsLibrary btkMathsLibrary btkDiffusionLibrary btkTractographyLibrary ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkODFParticleFilterTractography btkODFParticleFilterTractography.cxx
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
TARGET_LINK_LIBRARIES(btkODFParticleFilterTractography btkToolsLibrary ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkTensorParticleFilterTractography btkTensorParticleFilterTractography.cxx
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
TARGET_LINK_LIBRARIES(btkTensorParticleFilterTractography btkToolsLibrary ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkODFStreamlineTractography btkODFStreamlineTractography.cxx
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
TARGET_LINK_LIBRARIES(btkODFStreamlineTractography btkToolsLibrary ${ITK_LIBRARIES} vtkHybrid)

INSTALL(TARGETS
    btkODFParticleFilterTractography
    btkODFStreamlineTractography
    btkTensorParticleFilterTractography
    btkTractography
DESTINATION bin)

# ---- Segmentation ---------------------------------------------------------------------

ADD_EXECUTABLE(btkTissueSegmentation btkTissueSegmentation.cxx
    ${fbrain_SOURCE_DIR}/Code/Segmentation/btkTopologicalKMeans.h
)
TARGET_LINK_LIBRARIES(btkTissueSegmentation ${ITK_LIBRARIES})

INSTALL(TARGETS
    btkTissueSegmentation
DESTINATION bin)
