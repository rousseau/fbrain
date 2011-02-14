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

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Tractography)

ADD_EXECUTABLE(btkFiberODFStreamLineTractography btkFiberODFStreamLineTractography.cxx)
TARGET_LINK_LIBRARIES(btkFiberODFStreamLineTractography ITKIO vtkHybrid)

