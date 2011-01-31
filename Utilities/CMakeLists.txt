PROJECT(Utilities)
INCLUDE_REGULAR_EXPRESSION("^.*$")

ADD_EXECUTABLE(btkModifyImageUsingLookUpTable btkModifyImageUsingLookUpTable.cxx )
TARGET_LINK_LIBRARIES(btkModifyImageUsingLookUpTable ITKIO )

#ADD_EXECUTABLE(btkReorientImageToStandard btkReorientImageToStandard.cxx )
#TARGET_LINK_LIBRARIES(btkReorientImageToStandard ITKIO )

