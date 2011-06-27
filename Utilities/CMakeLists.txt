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
TARGET_LINK_LIBRARIES(btkModifyImageUsingLookUpTable ${ITK_LIBRARIES} )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkConvertGradientTable btkConvertGradientTable.cxx )
TARGET_LINK_LIBRARIES(btkConvertGradientTable ${ITK_LIBRARIES} )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkReorientImageToStandard btkReorientImageToStandard.cxx
        ${fbrain_SOURCE_DIR}/Code/InputOutput/btkLandmarksFileReader.cxx
)
TARGET_LINK_LIBRARIES(btkReorientImageToStandard ${ITK_LIBRARIES})
#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkReorientDiffusionSequenceToStandard btkReorientDiffusionSequenceToStandard.cxx
        ${fbrain_SOURCE_DIR}/Code/InputOutput/btkLandmarksFileReader.cxx
)

TARGET_LINK_LIBRARIES(btkReorientDiffusionSequenceToStandard ${ITK_LIBRARIES} ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkSetStandardCoorSystem btkSetStandardCoorSystem.cxx )
TARGET_LINK_LIBRARIES(btkSetStandardCoorSystem ${ITK_LIBRARIES} )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkPrintImageInfo btkPrintImageInfo.cxx )
TARGET_LINK_LIBRARIES(btkPrintImageInfo ${ITK_LIBRARIES} )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkRegisterDiffusionToAnatomicalData btkRegisterDiffusionToAnatomicalData.cxx )
TARGET_LINK_LIBRARIES(btkRegisterDiffusionToAnatomicalData ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNrrdToNifti btkNrrdToNifti.cxx
        ${fbrain_SOURCE_DIR}/Code/InputOutput/btkNrrdField.cxx)
TARGET_LINK_LIBRARIES(btkNrrdToNifti ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNiftiToNrrd btkNiftiToNrrd.cxx)
TARGET_LINK_LIBRARIES(btkNiftiToNrrd ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkExtractOneImageFromSequence btkExtractOneImageFromSequence.cxx)
TARGET_LINK_LIBRARIES(btkExtractOneImageFromSequence ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkCropImageUsingMask btkCropImageUsingMask.cxx)
TARGET_LINK_LIBRARIES(btkCropImageUsingMask ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkBinarizeLabels btkBinarizeLabels.cxx )
TARGET_LINK_LIBRARIES(btkBinarizeLabels ${ITK_LIBRARIES})

ADD_EXECUTABLE(btkResampleLabelsByInjection btkResampleLabelsByInjection.cxx )
TARGET_LINK_LIBRARIES(btkResampleLabelsByInjection ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkLowResolutionImageSimulator btkLowResolutionImageSimulator.cxx)
TARGET_LINK_LIBRARIES(btkLowResolutionImageSimulator ${ITK_LIBRARIES})

#ADD_EXECUTABLE(btkSimulateLowResolutionImage btkSimulateLowResolutionImage.cxx)
#TARGET_LINK_LIBRARIES(btkSimulateLowResolutionImage ${ITK_LIBRARIES})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkAverageImagesWithReference btkAverageImagesWithReference.cxx)
TARGET_LINK_LIBRARIES(btkAverageImagesWithReference ${ITK_LIBRARIES})

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
