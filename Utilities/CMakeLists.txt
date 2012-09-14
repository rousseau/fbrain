PROJECT(Utilities)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Reconstruction)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Diffusion)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/InputOutput)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Tools)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Transformations)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Segmentation)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Maths)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/ImageFilters)

ADD_EXECUTABLE(btkModifyImageUsingLookUpTable btkModifyImageUsingLookUpTable.cxx)
TARGET_LINK_LIBRARIES(btkModifyImageUsingLookUpTable ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkConvertGradientTable btkConvertGradientTable.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionGradientTable.h
 )
TARGET_LINK_LIBRARIES(btkConvertGradientTable ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkReorientImageToStandard btkReorientImageToStandard.cxx)
TARGET_LINK_LIBRARIES(btkReorientImageToStandard btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkReorientDiffusionSequenceToStandard btkReorientDiffusionSequenceToStandard.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionGradientTable.h
)
TARGET_LINK_LIBRARIES(btkReorientDiffusionSequenceToStandard btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSetStandardCoorSystem btkSetStandardCoorSystem.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionGradientTable.h
 )
TARGET_LINK_LIBRARIES(btkSetStandardCoorSystem ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkPrintImageInfo btkPrintImageInfo.cxx)
TARGET_LINK_LIBRARIES(btkPrintImageInfo ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkRegisterDiffusionToAnatomicalData btkRegisterDiffusionToAnatomicalData.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionGradientTable.h
 )
TARGET_LINK_LIBRARIES(btkRegisterDiffusionToAnatomicalData btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkNrrdToNifti btkNrrdToNifti.cxx)
TARGET_LINK_LIBRARIES(btkNrrdToNifti btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkNiftiToNrrd btkNiftiToNrrd.cxx)
TARGET_LINK_LIBRARIES(btkNiftiToNrrd btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkExtractOneImageFromSequence btkExtractOneImageFromSequence.cxx)
TARGET_LINK_LIBRARIES(btkExtractOneImageFromSequence ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkCropImageUsingMask btkCropImageUsingMask.cxx)
TARGET_LINK_LIBRARIES(btkCropImageUsingMask ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkBinarizeLabels btkBinarizeLabels.cxx )
TARGET_LINK_LIBRARIES(btkBinarizeLabels ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkResampleLabelsByInjection btkResampleLabelsByInjection.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkResampleLabelsByInjectionFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceTransform.h
 )
TARGET_LINK_LIBRARIES(btkResampleLabelsByInjection ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSimulateStandardViewFromIsotropicImage btkSimulateStandardViewFromIsotropicImage.cxx)
TARGET_LINK_LIBRARIES(btkSimulateStandardViewFromIsotropicImage ${ITK_LIBRARIES})

#ADD_EXECUTABLE(btkSimulateLowResolutionImage btkSimulateLowResolutionImage.cxx)
#TARGET_LINK_LIBRARIES(btkSimulateLowResolutionImage ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkAverageImagesWithReference btkAverageImagesWithReference.cxx)
TARGET_LINK_LIBRARIES(btkAverageImagesWithReference ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageInjection btkImageInjection.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkAffineRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRigidRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLowToHighImageResolutionMethod.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkResampleImageByInjectionFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkImageIntersectionCalculator.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceTransform.h
)
TARGET_LINK_LIBRARIES(btkImageInjection ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageSimilarity btkImageSimilarity.cxx)
TARGET_LINK_LIBRARIES(btkImageSimilarity ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkPSNR btkPSNR.cxx)
TARGET_LINK_LIBRARIES(btkPSNR ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkWeightedMean btkWeightedMean.cxx)
TARGET_LINK_LIBRARIES(btkWeightedMean ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSequenceNormalization btkSequenceNormalization.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionGradientTable.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkAffineRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkAffineRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Reconstruction/btkRigidRegistration.h
 )
TARGET_LINK_LIBRARIES(btkSequenceNormalization btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkRescaleIntensity btkRescaleIntensity.cxx)
TARGET_LINK_LIBRARIES(btkRescaleIntensity btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageHistogram btkImageHistogram.cxx)
TARGET_LINK_LIBRARIES(btkImageHistogram btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkReconstructionComparisonTool btkReconstructionComparisonTool.cxx)
TARGET_LINK_LIBRARIES(btkReconstructionComparisonTool ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkIteratedBackProjection btkIteratedBackProjection.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionRigidImageFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionAffineImageFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionManager.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceTransform.h
)
TARGET_LINK_LIBRARIES(btkIteratedBackProjection ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageResampling btkImageResampling.cxx)
TARGET_LINK_LIBRARIES(btkImageResampling ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkDifferentialBiasCorrection btkDifferentialBiasCorrection.cxx)
TARGET_LINK_LIBRARIES(btkDifferentialBiasCorrection ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkMajorityVoting btkMajorityVoting.cxx)
TARGET_LINK_LIBRARIES(btkMajorityVoting ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkProbabilityMapNormalization btkProbabilityMapNormalization.cxx)
TARGET_LINK_LIBRARIES(btkProbabilityMapNormalization ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageGaussianFilter btkImageGaussianFilter.cxx)
TARGET_LINK_LIBRARIES(btkImageGaussianFilter ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkBinarizeTissueProbabilityMaps btkBinarizeTissueProbabilityMaps.cxx)
TARGET_LINK_LIBRARIES(btkBinarizeTissueProbabilityMaps ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSimulateMotionSliceBySlice btkSimulateMotionSliceBySlice.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceRigidRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceAffineRegistration.h
)
TARGET_LINK_LIBRARIES(btkSimulateMotionSliceBySlice ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageSubtract btkImageSubtract.cxx)
TARGET_LINK_LIBRARIES(btkImageSubtract ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkApplyMaskToImage btkApplyMaskToImage.cxx)
TARGET_LINK_LIBRARIES(btkApplyMaskToImage ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkWarpTransformationToImage btkWarpTransformationToImage.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSliceBySliceTransform.h
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkWarpTransformToImageFilter.h
)
TARGET_LINK_LIBRARIES(btkWarpTransformationToImage ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkFCMClassification btkFCMClassification.cxx
    ${fbrain_SOURCE_DIR}/Code/Segmentation/btkFCMClassifier.h
)
TARGET_LINK_LIBRARIES(btkFCMClassification ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageMorphologicalClosing btkImageMorphologicalClosing.cxx)
TARGET_LINK_LIBRARIES(btkImageMorphologicalClosing ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageMorphologicalTopHat btkImageMorphologicalTopHat.cxx)
TARGET_LINK_LIBRARIES(btkImageMorphologicalTopHat ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkAverage3DImages btkAverage3DImages.cxx)
TARGET_LINK_LIBRARIES(btkAverage3DImages ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkComputeOverlap btkComputeOverlap.cxx)
TARGET_LINK_LIBRARIES(btkComputeOverlap ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkComputeChamferDistance btkComputeChamferDistance.cxx)
TARGET_LINK_LIBRARIES(btkComputeChamferDistance ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkMidwayHistogramEqualization btkMidwayHistogramEqualization.cxx
    ${fbrain_SOURCE_DIR}/Code/ImageFilters/btkMidwayImageEqualization.h
)
TARGET_LINK_LIBRARIES(btkMidwayHistogramEqualization ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkBinarizeMask btkBinarizeMask.cxx)
TARGET_LINK_LIBRARIES(btkBinarizeMask ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkExtractMaskUsingBoundingBox btkExtractMaskUsingBoundingBox.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkImageIntersectionCalculator.h
)
TARGET_LINK_LIBRARIES(btkExtractMaskUsingBoundingBox ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkTranslateImageOverTemplate btkTranslateImageOverTemplate.cxx)
TARGET_LINK_LIBRARIES(btkTranslateImageOverTemplate ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkExtendMaskToImage btkExtendMaskToImage.cxx)
TARGET_LINK_LIBRARIES(btkExtendMaskToImage ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkWeightedSumOfAffineTransforms btkWeightedSumOfAffineTransforms.cxx)
TARGET_LINK_LIBRARIES(btkWeightedSumOfAffineTransforms btkToolsLibrary btkMathsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkInverseDisplacementField btkInverseDisplacementField.cxx)
TARGET_LINK_LIBRARIES(btkInverseDisplacementField ${ITK_LIBRARIES})

INSTALL(TARGETS btkModifyImageUsingLookUpTable
                btkConvertGradientTable
                btkReorientImageToStandard
                btkReorientDiffusionSequenceToStandard
                btkSetStandardCoorSystem
                btkPrintImageInfo
                btkRegisterDiffusionToAnatomicalData
                btkNrrdToNifti
                btkNiftiToNrrd
                btkExtractOneImageFromSequence
                btkCropImageUsingMask
                btkBinarizeLabels
                btkResampleLabelsByInjection
                btkSimulateStandardViewFromIsotropicImage
                btkAverageImagesWithReference
                btkImageInjection
                btkImageSimilarity
                btkPSNR
                btkWeightedMean
                btkSequenceNormalization
                btkReconstructionComparisonTool
                btkIteratedBackProjection
                btkImageResampling
                btkDifferentialBiasCorrection
                btkMajorityVoting
                btkProbabilityMapNormalization
                btkImageGaussianFilter
                btkBinarizeTissueProbabilityMaps
                btkSimulateMotionSliceBySlice
                btkImageSubtract
                btkApplyMaskToImage
                btkWarpTransformationToImage
                btkFCMClassification
                btkImageMorphologicalClosing
                btkImageMorphologicalTopHat
                btkAverage3DImages
                btkComputeOverlap
                btkComputeChamferDistance
                btkMidwayHistogramEqualization
                btkBinarizeMask
                btkExtractMaskUsingBoundingBox
                btkTranslateImageOverTemplate
                btkExtendMaskToImage
                btkImageHistogram
                btkRescaleIntensity
                btkWeightedSumOfAffineTransforms
                btkInverseDisplacementField
DESTINATION bin)
