PROJECT(Utilities)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})
INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/DTI)

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkModifyImageUsingLookUpTable btkModifyImageUsingLookUpTable.cxx )
TARGET_LINK_LIBRARIES(btkModifyImageUsingLookUpTable ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkReorientImageToStandard btkReorientImageToStandard.cxx )
TARGET_LINK_LIBRARIES(btkReorientImageToStandard ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkReorientDiffusionSequenceToStandard btkReorientDiffusionSequenceToStandard.cxx )
TARGET_LINK_LIBRARIES(btkReorientDiffusionSequenceToStandard ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkSetStandardCoorSystem btkSetStandardCoorSystem.cxx )
TARGET_LINK_LIBRARIES(btkSetStandardCoorSystem ITKIO )

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkNrrdToNifti btkNrrdToNifti.cxx)
TARGET_LINK_LIBRARIES(btkNrrdToNifti ITKIO)
