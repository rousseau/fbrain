PROJECT(Utilities)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/DTI)
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Various)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkModifyImageUsingLookUpTable btkModifyImageUsingLookUpTable.cxx )
TARGET_LINK_LIBRARIES(btkModifyImageUsingLookUpTable ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkConvertGradientTable btkConvertGradientTable.cxx )
TARGET_LINK_LIBRARIES(btkConvertGradientTable ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkReorientImageToStandard btkReorientImageToStandard.cxx
        ${fbrain_SOURCE_DIR}/Code/Various/btkLandmarksFileReader.cxx
)
TARGET_LINK_LIBRARIES(btkReorientImageToStandard ITKCommon ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkReorientDiffusionSequenceToStandard btkReorientDiffusionSequenceToStandard.cxx
        ${fbrain_SOURCE_DIR}/Code/Various/btkLandmarksFileReader.cxx
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

ADD_EXECUTABLE(btkNrrdToNifti btkNrrdToNifti.cxx)
TARGET_LINK_LIBRARIES(btkNrrdToNifti ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNiftiToNrrd btkNiftiToNrrd.cxx)
TARGET_LINK_LIBRARIES(btkNiftiToNrrd ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkDWIImageExtractor btkDWIImageExtractor.cxx)
TARGET_LINK_LIBRARIES(btkDWIImageExtractor ITKIO)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkCropImageUsingMask btkCropImageUsingMask.cxx)
TARGET_LINK_LIBRARIES(btkCropImageUsingMask ITKIO)

#------------------------------------------------------------------------------

INSTALL(TARGETS btkModifyImageUsingLookUpTable
                btkReorientImageToStandard
                btkReorientDiffusionSequenceToStandard
                btkSetStandardCoorSystem
                btkNiftiToNrrd
                btkNrrdToNifti
                btkDWIImageExtractor
		btkCropImageUsingMask
    btkCropImageUsingMask
        DESTINATION bin)
