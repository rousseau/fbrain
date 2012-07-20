PROJECT(Code)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#---- Directories -----------------------------------------------------------------------

SET(TOOLS_LIBRARY_SOURCE_DIR ${Code_SOURCE_DIR}/Tools)
SET(MATHS_LIBRARY_SOURCE_DIR ${Code_SOURCE_DIR}/Maths)

IF(BUILD_TESTING)
    SET(MATHS_TESTS_SOURCE_DIR ${MATHS_LIBRARY_SOURCE_DIR}/Tests)
ENDIF(BUILD_TESTING)

INCLUDE_DIRECTORIES(${MATHS_LIBRARY_SOURCE_DIR} ${TOOLS_LIBRARY_SOURCE_DIR})

#---- Maths library ---------------------------------------------------------------------

ADD_LIBRARY(btkMathsLibrary STATIC
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkMacro.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkMathFunctions.cxx ${MATHS_LIBRARY_SOURCE_DIR}/btkMathFunctions.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalDirection.cxx ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalDirection.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkLegendrePolynomial.cxx ${MATHS_LIBRARY_SOURCE_DIR}/btkLegendrePolynomial.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalHarmonics.cxx ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalHarmonics.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkMatrixOperations.cxx ${MATHS_LIBRARY_SOURCE_DIR}/btkMatrixOperations.h
)

IF(BUILD_TESTING)
    ADD_EXECUTABLE(btkMathsLibraryTestsApp ${MATHS_TESTS_SOURCE_DIR}/btkMathsLibraryTestsApp.cxx)
    TARGET_LINK_LIBRARIES(btkMathsLibraryTestsApp btkMathsLibrary ${CPPUNIT_LIBRARY})
    ADD_TEST(btkMathsLibraryTests btkMathsLibraryTestsApp)
ENDIF(BUILD_TESTING)
