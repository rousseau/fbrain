PROJECT(Tests)
INCLUDE_REGULAR_EXPRESSION("^.*$")


ADD_EXECUTABLE(btkDiffusionSequenceTestApp ${Tests_SOURCE_DIR}/btkDiffusionSequenceTest.cxx)
ADD_TEST(btkDiffusionSequenceTest ${Tests_BINARY_DIR}/btkDiffusionSequenceTestApp)
