PROJECT(Utilities)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Reconstruction)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/DTI)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/InputOutput)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Reconstruction)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkModifyImageUsingLookUpTable btkModifyImageUsingLookUpTable.cxx )
TARGET_LINK_LIBRARIES(btkModifyImageUsingLookUpTable ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkConvertGradientTable btkConvertGradientTable.cxx )
TARGET_LINK_LIBRARIES(btkConvertGradientTable ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkReorientImageToStandard btkReorientImageToStandard.cxx
        ${fbrain_SOURCE_DIR}/Code/InputOutput/btkLandmarksFileReader.cxx
)
TARGET_LINK_LIBRARIES(btkReorientImageToStandard ITKCommon ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkReorientDiffusionSequenceToStandard btkReorientDiffusionSequenceToStandard.cxx
        ${fbrain_SOURCE_DIR}/Code/InputOutput/btkLandmarksFileReader.cxx
)

TARGET_LINK_LIBRARIES(btkReorientDiffusionSequenceToStandard ITKCommon ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkSetStandardCoorSystem btkSetStandardCoorSystem.cxx )
TARGET_LINK_LIBRARIES(btkSetStandardCoorSystem ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkPrintImageInfo btkPrintImageInfo.cxx )
TARGET_LINK_LIBRARIES(btkPrintImageInfo ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkRegisterDiffusionToAnatomicalData btkRegisterDiffusionToAnatomicalData.cxx )
TARGET_LINK_LIBRARIES(btkRegisterDiffusionToAnatomicalData ITKCommon ITKIO ITKNumerics)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNrrdToNifti btkNrrdToNifti.cxx
        ${fbrain_SOURCE_DIR}/Code/InputOutput/btkNrrdField.cxx)
TARGET_LINK_LIBRARIES(btkNrrdToNifti ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNiftiToNrrd btkNiftiToNrrd.cxx)
TARGET_LINK_LIBRARIES(btkNiftiToNrrd ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkExtractOneImageFromSequence btkExtractOneImageFromSequence.cxx)
TARGET_LINK_LIBRARIES(btkExtractOneImageFromSequence ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkCropImageUsingMask btkCropImageUsingMask.cxx)
TARGET_LINK_LIBRARIES(btkCropImageUsingMask ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkBinarizeLabels btkBinarizeLabels.cxx )
TARGET_LINK_LIBRARIES(btkBinarizeLabels ITKIO ITKNumerics ITKStatistics)

ADD_EXECUTABLE(btkResampleLabelsByInjection btkResampleLabelsByInjection.cxx )
TARGET_LINK_LIBRARIES(btkResampleLabelsByInjection ITKIO ITKNumerics )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkLowResolutionImageSimulator btkLowResolutionImageSimulator.cxx)
TARGET_LINK_LIBRARIES(btkLowResolutionImageSimulator ITKIO)

ADD_EXECUTABLE(btkSimulateLowResolutionImage btkSimulateLowResolutionImage.cxx)
TARGET_LINK_LIBRARIES(btkSimulateLowResolutionImage ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkAverageImagesWithReference btkAverageImagesWithReference.cxx)
TARGET_LINK_LIBRARIES(btkAverageImagesWithReference ITKIO)

#------------------------------------------------------------------------------

INSTALL(TARGETS btkModifyImageUsingLookUpTable
                btkReorientImageToStandard
                btkReorientDiffusionSequenceToStandard
                btkSetStandardCoorSystem
                btkNiftiToNrrd
                btkNrrdToNifti
                btkExtractOneImageFromSequence
                btkCropImageUsingMask
        DESTINATION bin)
