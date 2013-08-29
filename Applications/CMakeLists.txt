PROJECT(Applications)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INCLUDE_DIRECTORIES(
    ${TCLAP_DIRECTORY}
    ${Code_SOURCE_DIR}/Tools
    ${Code_SOURCE_DIR}/Maths
    ${Code_SOURCE_DIR}/Transformations
    ${Code_SOURCE_DIR}/Visualization

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
    ${Code_SOURCE_DIR}/Maths
    ${Code_SOURCE_DIR}/Visualization
    ${Code_SOURCE_DIR}/ImageFilters
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
    ${fbrain_SOURCE_DIR}/Code/Maths/btkPSF.h
    ${fbrain_SOURCE_DIR}/Code/Maths/btkPSF.cxx
    ${fbrain_SOURCE_DIR}/Code/Maths/btkGaussianPSF.h
    ${fbrain_SOURCE_DIR}/Code/Maths/btkGaussianPSF.cxx
)
TARGET_LINK_LIBRARIES(btkSuperResolution ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSuperResolutionV2 btkSuperResolutionV2.cxx
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
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.h
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkMseCostFunction.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceRigidRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceAffineRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Maths/btkPSF.h
    ${fbrain_SOURCE_DIR}/Code/Maths/btkPSF.cxx
    ${fbrain_SOURCE_DIR}/Code/Maths/btkGaussianPSF.h
    ${fbrain_SOURCE_DIR}/Code/Maths/btkGaussianPSF.cxx
    ${fbrain_SOURCE_DIR}/Code/Maths/btkBoxCarPSF.h
    ${fbrain_SOURCE_DIR}/Code/Maths/btkBoxCarPSF.cxx
    ${fbrain_SOURCE_DIR}/Code/Maths/btkSincPSF.h
    ${fbrain_SOURCE_DIR}/Code/Maths/btkSincPSF.cxx
)
TARGET_LINK_LIBRARIES(btkSuperResolutionV2 btkToolsLibrary btkMathsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkGroupwiseS2SDistortionCorrection btkGroupwiseS2SDistortionCorrection.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolation.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkNumerical.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkGroupwiseS2SDistortionCorrection.h
)
TARGET_LINK_LIBRARIES(btkGroupwiseS2SDistortionCorrection btkToolsLibrary ${ITK_LIBRARIES} ${ANN_LIBRARY})


#ADD_EXECUTABLE(btk3DSuperResolution btk3DSuperResolution.cxx
 #   ${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkSuperResolutionFilter.cxx
  #  ${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkSuperResolutionFilter.h
   # ${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkBiasCorrectionFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkBiasCorrectionFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLowToHighImageResolutionMethod.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionReconstructionFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionIBPFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionSRFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionReconstructionFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionIBPFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkHighResolutionSRFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionSliceBySliceAffineFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionSliceBySliceAffineFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionSliceBySliceEulerFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrectionSliceBySliceEulerFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrection3DAffineFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrection3DAffineFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrection3DEulerFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/SuperResolution/btkMotionCorrection3DEulerFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSimulateLRImageFilter.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.cxx
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkCreateHRMaskFilter.h
    #${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    #${fbrain_SOURCE_DIR}/Code/Registration/btkMseCostFunction.h
    #${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceRigidRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceAffineRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkMseCostFunction.h
#)

#TARGET_LINK_LIBRARIES(btk3DSuperResolution btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkRBFInterpolationS2S btkRBFInterpolationS2S.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolation.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkNumerical.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolateImageFunctionS2S.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionGradientTable.h
)

TARGET_LINK_LIBRARIES(btkRBFInterpolationS2S btkToolsLibrary  btkDiffusionLibrary ${ITK_LIBRARIES} ${ANN_LIBRARY})

#ADD_EXECUTABLE(btkIntersectionRegistration btkIntersectionRegistration.cxx
#${fbrain_SOURCE_DIR}/Code/Tools/btkImageHelper.h
#${fbrain_SOURCE_DIR}/Code/Registration/btkSlicesIntersectionVNLCostFunction.hxx
#${fbrain_SOURCE_DIR}/Code/Registration/btkSlicesIntersectionITKCostFunction.hxx
#${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLowToHighResolutionFilter.hxx
#${fbrain_SOURCE_DIR}/Code/Registration/btkOptimizer.cxx ${fbrain_SOURCE_DIR}/Code/Registration/btkOptimizer.h
#${fbrain_SOURCE_DIR}/Code/Registration/btkSimulatedAnnealingOptimizer.cxx ${fbrain_SOURCE_DIR}/Code/Registration/btkSimulatedAnnealingOptimizer.h
#${fbrain_SOURCE_DIR}/Code/Registration/btkSmartStepGradientDescentOptimizer.cxx ${fbrain_SOURCE_DIR}/Code/Registration/btkSmartStepGradientDescentOptimizer.h
#${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSRHMatrixComputation.hxx
#${fbrain_SOURCE_DIR}/Code/Maths/btkPSF.h
#${fbrain_SOURCE_DIR}/Code/Maths/btkPSF.cxx
#${fbrain_SOURCE_DIR}/Code/Maths/btkGaussianPSF.h
#${fbrain_SOURCE_DIR}/Code/Maths/btkGaussianPSF.cxx
#${fbrain_SOURCE_DIR}/Code/Maths/btkBoxCarPSF.h
#${fbrain_SOURCE_DIR}/Code/Maths/btkBoxCarPSF.cxx)

#TARGET_LINK_LIBRARIES(btkIntersectionRegistration btkToolsLibrary btkMathsLibrary ${ITK_LIBRARIES} vtkHybrid)

INSTALL(TARGETS btkImageReconstruction
                btkSuperResolution
                btkSuperResolutionV2
                btkGroupwiseS2SDistortionCorrection
                btkRBFInterpolationS2S
               #btkIntersectionRegistration
        DESTINATION bin)



# ---- Tractography ---------------------------------------------------------------------

INCLUDE_DIRECTORIES(
    ${Code_SOURCE_DIR}/Denoising
    ${Code_SOURCE_DIR}/Diffusion
    ${Code_SOURCE_DIR}/Tractography
)

ADD_EXECUTABLE(btkTractography btkTractography.cxx
    ${Code_SOURCE_DIR}/Denoising/btkPseudoResidualsVarianceCalculator.h
)
TARGET_LINK_LIBRARIES(btkTractography btkToolsLibrary btkMathsLibrary btkDiffusionLibrary btkTractographyLibrary btkDenoisingLibrary ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkAtlasBasedFiberLabeling btkAtlasBasedFiberLabeling.cxx
)
TARGET_LINK_LIBRARIES(btkAtlasBasedFiberLabeling ${ITK_LIBRARIES} vtkHybrid)

INSTALL(TARGETS
    btkTractography
    btkAtlasBasedFiberLabeling
DESTINATION bin)

# ---- Segmentation ---------------------------------------------------------------------

ADD_EXECUTABLE(btkTissueSegmentation btkTissueSegmentation.cxx
    ${fbrain_SOURCE_DIR}/Code/Segmentation/btkTopologicalKMeans.h
)
TARGET_LINK_LIBRARIES(btkTissueSegmentation ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkLabelPropagation btkLabelPropagation.cxx
    ${fbrain_SOURCE_DIR}/Code/Segmentation/btkLabelPropagationTool.h
)
TARGET_LINK_LIBRARIES(btkLabelPropagation ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSSL btkSSL.cxx)
TARGET_LINK_LIBRARIES(btkSSL btkToolsLibrary ${ITK_LIBRARIES})

INSTALL(TARGETS
    btkTissueSegmentation
    btkLabelPropagation
    btkSSL
DESTINATION bin)

# ---- Clustering -----------------------------------------------------------------------

ADD_EXECUTABLE( btkProbabilisticSegmentationMapBasedClustering  btkProbabilisticSegmentationMapBasedClustering.cxx)
TARGET_LINK_LIBRARIES(btkProbabilisticSegmentationMapBasedClustering ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE( btkDistanceFibersKMeansApproximationClustering  btkDistanceFibersKMeansApproximationClustering.cxx)
TARGET_LINK_LIBRARIES(btkDistanceFibersKMeansApproximationClustering ${ITK_LIBRARIES} vtkHybrid vtkInfovis)

INSTALL(TARGETS 
    btkProbabilisticSegmentationMapBasedClustering
    btkDistanceFibersKMeansApproximationClustering
DESTINATION bin)

# ---- Outliers correction -----------------------------------------------------------------------
ADD_EXECUTABLE(btkOutliersCorrection btkOutliersCorrection.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLUDecomposition.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRBFInterpolateImageFunctionS2S.h
)
TARGET_LINK_LIBRARIES(btkOutliersCorrection btkToolsLibrary btkDiffusionLibrary ${ITK_LIBRARIES} ${ANN_LIBRARY})

INSTALL(TARGETS
    btkOutliersCorrection
DESTINATION bin)


# ---- Feature selection ----------------------------------------------------------------

include_directories(
    ${Code_SOURCE_DIR}/FeatureSelection
)

add_executable(btkFeatureSelection btkFeatureSelection.cxx)
target_link_libraries(btkFeatureSelection ${ITK_LIBRARIES} btkFeatureSelectionLibrary)

install(TARGETS
    btkFeatureSelection
DESTINATION bin)

# ---- Dwi registrtion correction -----------------------------------------------------------------------

ADD_EXECUTABLE(btkDwiReconstruction btkDwiReconstruction.cxx
${Code_SOURCE_DIR}/Reconstruction/btkDwiReconstructionFilter.h
${Code_SOURCE_DIR}/Reconstruction/btkDwiReconstructionFilter.cxx)

TARGET_LINK_LIBRARIES(btkDwiReconstruction btkToolsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})


INSTALL(TARGETS
    btkDwiReconstruction
DESTINATION bin)
