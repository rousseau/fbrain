PROJECT(Utilities)
INCLUDE_REGULAR_EXPRESSION("^.*$")

ADD_EXECUTABLE(ModifyImageUsingLookUpTable ModifyImageUsingLookUpTable.cxx )
TARGET_LINK_LIBRARIES(ModifyImageUsingLookUpTable  ITKIO boost_program_options)


