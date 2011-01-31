PROJECT(Utilities)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${TCLAP_DIRECTORY})

#------------------------------------------------------------------------------

ADD_EXECUTABLE(btkModifyImageUsingLookUpTable btkModifyImageUsingLookUpTable.cxx )
TARGET_LINK_LIBRARIES(btkModifyImageUsingLookUpTable ITKIO )


