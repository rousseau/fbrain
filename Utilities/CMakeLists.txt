PROJECT(Utilities)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})
INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Reconstruction)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Diffusion)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Tools)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Transformations)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Segmentation)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Maths)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/ImageFilters)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Registration)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Tractography)
INCLUDE_DIRECTORIES(${Code_SOURCE_DIR}/Visualization)
include_directories(${Code_SOURCE_DIR}/FeatureSelection)

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
TARGET_LINK_LIBRARIES(btkPrintImageInfo ${ITK_LIBRARIES} btkToolsLibrary btkDiffusionLibrary)

ADD_EXECUTABLE(btkRegisterDiffusionToAnatomicalData btkRegisterDiffusionToAnatomicalData.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionGradientTable.h
 )
TARGET_LINK_LIBRARIES(btkRegisterDiffusionToAnatomicalData btkToolsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkNrrdToNifti btkNrrdToNifti.cxx)
TARGET_LINK_LIBRARIES(btkNrrdToNifti btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkNiftiToNrrd btkNiftiToNrrd.cxx)
TARGET_LINK_LIBRARIES(btkNiftiToNrrd btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkExtractOneImageFromSequence btkExtractOneImageFromSequence.cxx)
TARGET_LINK_LIBRARIES(btkExtractOneImageFromSequence ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkCropImageUsingMask btkCropImageUsingMask.cxx)
TARGET_LINK_LIBRARIES(btkCropImageUsingMask btkToolsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkBinarizeLabels btkBinarizeLabels.cxx )
TARGET_LINK_LIBRARIES(btkBinarizeLabels ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkResampleLabelsByInjection btkResampleLabelsByInjection.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkResampleLabelsByInjectionFilter.h
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
 )
TARGET_LINK_LIBRARIES(btkResampleLabelsByInjection ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSimulateStandardViewFromIsotropicImage btkSimulateStandardViewFromIsotropicImage.cxx)
TARGET_LINK_LIBRARIES(btkSimulateStandardViewFromIsotropicImage ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSimulateLowResolutionImage btkSimulateLowResolutionImage.cxx
${fbrain_SOURCE_DIR}/Code/Transformations/btkRandomSliceBySliceTransformGenerator.h
${fbrain_SOURCE_DIR}/Code/Transformations/btkRandomSliceBySliceTransformGenerator.cxx)
TARGET_LINK_LIBRARIES(btkSimulateLowResolutionImage btkMathsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkAverageImagesWithReference btkAverageImagesWithReference.cxx)
TARGET_LINK_LIBRARIES(btkAverageImagesWithReference ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageInjection btkImageInjection.cxx
    ${fbrain_SOURCE_DIR}/Code/Registration/btkRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkAffineRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkRigidRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkLowToHighImageResolutionMethod.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkResampleImageByInjectionFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkImageIntersectionCalculator.h
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
)
TARGET_LINK_LIBRARIES(btkImageInjection ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageSimilarity btkImageSimilarity.cxx)
TARGET_LINK_LIBRARIES(btkImageSimilarity ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkPSNR btkPSNR.cxx)
TARGET_LINK_LIBRARIES(btkPSNR ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkWeightedSumOfImages btkWeightedSumOfImages.cxx
    ${fbrain_SOURCE_DIR}/Code/ImageFilters/btkResampleImagesToBiggestImageFilter.h
    ${fbrain_SOURCE_DIR}/Code/ImageFilters/btkWeightedSumOfImagesFilter.h)
TARGET_LINK_LIBRARIES(btkWeightedSumOfImages ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSequenceNormalization btkSequenceNormalization.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionGradientTable.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkAffineRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Registration/btkRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Registration/btkAffineRegistration.h
    #${fbrain_SOURCE_DIR}/Code/Registration/btkRigidRegistration.h
 )
TARGET_LINK_LIBRARIES(btkSequenceNormalization btkToolsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkRescaleIntensity btkRescaleIntensity.cxx)
TARGET_LINK_LIBRARIES(btkRescaleIntensity btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageHistogram btkImageHistogram.cxx)
TARGET_LINK_LIBRARIES(btkImageHistogram btkToolsLibrary btkMathsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkReconstructionComparisonTool btkReconstructionComparisonTool.cxx)
TARGET_LINK_LIBRARIES(btkReconstructionComparisonTool ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkIteratedBackProjection btkIteratedBackProjection.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionRigidImageFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionAffineImageFilter.h
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkSuperResolutionManager.h
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
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
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceRigidRegistration.h
    ${fbrain_SOURCE_DIR}/Code/Registration/btkSliceBySliceAffineRegistration.h
)
TARGET_LINK_LIBRARIES(btkSimulateMotionSliceBySlice btkMathsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkImageSubtract btkImageSubtract.cxx)
TARGET_LINK_LIBRARIES(btkImageSubtract ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkApplyMaskToImage btkApplyMaskToImage.cxx)
TARGET_LINK_LIBRARIES(btkApplyMaskToImage ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkApplyTransformationToImage btkApplyTransformationToImage.cxx
    ${fbrain_SOURCE_DIR}/Code/Transformations/btkSliceBySliceTransform.h
    ${fbrain_SOURCE_DIR}/Code/ImageFilters/btkApplyTransformToImageFilter.h
)
TARGET_LINK_LIBRARIES(btkApplyTransformationToImage ${ITK_LIBRARIES})

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
TARGET_LINK_LIBRARIES(btkMidwayHistogramEqualization btkMathsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkHistogramMatching btkHistogramMatching.cxx)
TARGET_LINK_LIBRARIES(btkHistogramMatching ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkDistanceBetweenBinaryImages btkDistanceBetweenBinaryImages.cxx)
TARGET_LINK_LIBRARIES(btkDistanceBetweenBinaryImages ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkBinarizeMask btkBinarizeMask.cxx)
TARGET_LINK_LIBRARIES(btkBinarizeMask ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkExtractMaskUsingBoundingBox btkExtractMaskUsingBoundingBox.cxx
    ${fbrain_SOURCE_DIR}/Code/Reconstruction/btkImageIntersectionCalculator.h
)
TARGET_LINK_LIBRARIES(btkExtractMaskUsingBoundingBox ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkTranslateImageOverTemplate btkTranslateImageOverTemplate.cxx)
TARGET_LINK_LIBRARIES(btkTranslateImageOverTemplate ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkWeightedSumOfAffineTransforms btkWeightedSumOfAffineTransforms.cxx)
TARGET_LINK_LIBRARIES(btkWeightedSumOfAffineTransforms btkToolsLibrary btkMathsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkInverseDisplacementField btkInverseDisplacementField.cxx)
TARGET_LINK_LIBRARIES(btkInverseDisplacementField ${ITK_LIBRARIES})

ADD_EXECUTABLE( btkTractDensityMap  btkTractDensityMap.cxx)
TARGET_LINK_LIBRARIES(btkTractDensityMap ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkColorFiberTractsByOrientation btkColorFiberTractsByOrientation.cxx)
TARGET_LINK_LIBRARIES(btkColorFiberTractsByOrientation btkTractographyLibrary vtkHybrid)

ADD_EXECUTABLE(btkComputeSoftMaskUsingOrthogonalImages btkComputeSoftMaskUsingOrthogonalImages.cxx)
TARGET_LINK_LIBRARIES(btkComputeSoftMaskUsingOrthogonalImages ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkLocalHistogramMatching btkLocalHistogramMatching.cxx)
TARGET_LINK_LIBRARIES(btkLocalHistogramMatching btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkWeightImageOfPatches btkWeightImageOfPatches.cxx)
TARGET_LINK_LIBRARIES(btkWeightImageOfPatches btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE( btkTriangleMesh btkTriangleMesh.cxx )
TARGET_LINK_LIBRARIES(btkTriangleMesh vtkHybrid ${VTK_LIBRARIES})

ADD_EXECUTABLE( btkDiffusionScalarMeasurement btkDiffusionScalarMeasurement.cxx)
TARGET_LINK_LIBRARIES(btkDiffusionScalarMeasurement btkToolsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkS2SSimilarity btkS2SSimilarity.cxx)
TARGET_LINK_LIBRARIES(btkS2SSimilarity btkToolsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkComputeCurvature btkComputeCurvature.cxx)
TARGET_LINK_LIBRARIES(btkComputeCurvature btkMathsLibrary ${VTK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkConvertToTXT btkConvertToTXT.cxx)
TARGET_LINK_LIBRARIES(btkConvertToTXT ${VTK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkApplyRegionGrow btkApplyRegionGrow.cxx)
TARGET_LINK_LIBRARIES(btkApplyRegionGrow btkToolsLibrary btkMathsLibrary ${VTK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkProjectPointsOnPolydata btkProjectPointsOnPolydata.cxx)
TARGET_LINK_LIBRARIES(btkProjectPointsOnPolydata ${VTK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkSelectSeedPoint btkSelectSeedPoint.cxx)
TARGET_LINK_LIBRARIES(btkSelectSeedPoint ${VTK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkLikelihoodOfSequenceEstimation btkLikelihoodOfSequenceEstimation.cxx)
TARGET_LINK_LIBRARIES(btkLikelihoodOfSequenceEstimation ${ITK_LIBRARIES})

add_executable(btkComputeFeatureSelectionResiduals btkComputeFeatureSelectionResiduals.cxx)
target_link_libraries(btkComputeFeatureSelectionResiduals btkFeatureSelectionLibrary ${ITK_LIBRARIES})

add_executable(btkComputeCoefficientOfDeterminationForAtlas btkComputeCoefficientOfDeterminationForAtlas.cxx)
target_link_libraries(btkComputeCoefficientOfDeterminationForAtlas btkToolsLibrary ${ITK_LIBRARIES})


ADD_EXECUTABLE( btkGenerateVtkFileFromFiberDataTextFiles btkGenerateVtkFileFromFiberDataTextFiles.cxx)
TARGET_LINK_LIBRARIES(btkGenerateVtkFileFromFiberDataTextFiles ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE( btkGenerateSimulatedFiber  btkGenerateSimulatedFiber.cxx)
TARGET_LINK_LIBRARIES(btkGenerateSimulatedFiber ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE( btkCreate4DEmptyImageFromReferenceImage btkCreate4DEmptyImageFromReferenceImage.cxx)
TARGET_LINK_LIBRARIES(btkCreate4DEmptyImageFromReferenceImage ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE( btkQuantifyClusteredRegions btkQuantifyClusteredRegions.cxx)
TARGET_LINK_LIBRARIES(btkQuantifyClusteredRegions ${ITK_LIBRARIES} vtkHybrid)

add_executable(btkApplyTransformToDiffusionSequence btkApplyTransformToDiffusionSequence.cxx)
target_link_libraries(btkApplyTransformToDiffusionSequence btkToolsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})


add_executable(btkApplyAffineTransformToFibers btkApplyAffineTransformToFibers.cxx)
target_link_libraries(btkApplyAffineTransformToFibers btkToolsLibrary ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkResolutionSampler btkResolutionSampler.cxx)
TARGET_LINK_LIBRARIES(btkResolutionSampler  ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkSTAPLE btkSTAPLE.cxx)
TARGET_LINK_LIBRARIES(btkSTAPLE  ${ITK_LIBRARIES})

add_executable(btkCreateMarkupListFromSelectedFeatures btkCreateMarkupListFromSelectedFeatures.cxx)
target_link_libraries(btkCreateMarkupListFromSelectedFeatures btkFeatureSelectionLibrary ${ITK_LIBRARIES})

add_executable(btkFeaturesReconstruction btkFeaturesReconstruction.cxx)
target_link_libraries(btkFeaturesReconstruction btkFeatureSelectionLibrary ${ITK_LIBRARIES})

add_executable(btkFeaturesStatistics btkFeaturesStatistics.cxx)
target_link_libraries(btkFeaturesStatistics btkFeatureSelectionLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkComputeStatisticsOverMultipleROI btkComputeStatisticsOverMultipleROI.cxx)
TARGET_LINK_LIBRARIES(btkComputeStatisticsOverMultipleROI  btkToolsLibrary ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkComputeJacobian btkComputeJacobian.cxx)
TARGET_LINK_LIBRARIES(btkComputeJacobian  ${ITK_LIBRARIES})

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
                btkSimulateLowResolutionImage
                btkAverageImagesWithReference
                btkImageInjection
                btkImageSimilarity
                btkPSNR
                btkWeightedSumOfImages
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
                btkApplyTransformationToImage
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
                btkImageHistogram
                btkRescaleIntensity
                btkHistogramMatching
                btkWeightedSumOfAffineTransforms
                btkInverseDisplacementField
                btkDistanceBetweenBinaryImages
                btkTractDensityMap
                btkColorFiberTractsByOrientation
                btkComputeSoftMaskUsingOrthogonalImages
                btkLocalHistogramMatching
                btkWeightImageOfPatches
                btkTriangleMesh
                btkDiffusionScalarMeasurement
                btkS2SSimilarity
                btkComputeCurvature
                btkConvertToTXT
                btkApplyRegionGrow
                btkProjectPointsOnPolydata
                btkSelectSeedPoint
                btkLikelihoodOfSequenceEstimation
                btkComputeFeatureSelectionResiduals
                btkComputeCoefficientOfDeterminationForAtlas
                btkGenerateVtkFileFromFiberDataTextFiles
                btkGenerateSimulatedFiber
				btkCreate4DEmptyImageFromReferenceImage
				btkQuantifyClusteredRegions
                btkApplyTransformToDiffusionSequence
                btkApplyAffineTransformToFibers
                btkResolutionSampler
                btkSTAPLE
                btkCreateMarkupListFromSelectedFeatures
                btkFeaturesReconstruction
                btkFeaturesStatistics
                btkComputeStatisticsOverMultipleROI
                btkComputeJacobian
        DESTINATION bin)

