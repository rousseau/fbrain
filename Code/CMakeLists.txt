PROJECT(Code)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#---- Directories -----------------------------------------------------------------------

SET(TOOLS_LIBRARY_SOURCE_DIR     ${Code_SOURCE_DIR}/Tools)
SET(MATHS_LIBRARY_SOURCE_DIR     ${Code_SOURCE_DIR}/Maths)
SET(DIFFUSION_LIBRARY_SOURCE_DIR ${Code_SOURCE_DIR}/Diffusion)

IF(BUILD_TESTING)
    SET(MATHS_TESTS_SOURCE_DIR ${MATHS_LIBRARY_SOURCE_DIR}/Tests)
ENDIF(BUILD_TESTING)

INCLUDE_DIRECTORIES(
    ${MATHS_LIBRARY_SOURCE_DIR}
    ${TOOLS_LIBRARY_SOURCE_DIR}
    ${DIFFUSION_LIBRARY_SOURCE_DIR}
)

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
    ADD_EXECUTABLE(btkMathsLibraryTestsApp ${MATHS_TESTS_SOURCE_DIR}/btkMathsLibraryTestsApp.cxx
        ${MATHS_TESTS_SOURCE_DIR}/btkMatrixOperationsTest.cxx ${MATHS_TESTS_SOURCE_DIR}/btkMatrixOperationsTest.h
        ${MATHS_TESTS_SOURCE_DIR}/btkLegendrePolynomialTest.cxx ${MATHS_TESTS_SOURCE_DIR}/btkLegendrePolynomialTest.h
        ${MATHS_TESTS_SOURCE_DIR}/btkMathFunctionsTest.cxx ${MATHS_TESTS_SOURCE_DIR}/btkMathFunctionsTest.h
        ${MATHS_TESTS_SOURCE_DIR}/btkSphericalHarmonicsTest.cxx ${MATHS_TESTS_SOURCE_DIR}/btkSphericalHarmonicsTest.h
    )
    TARGET_LINK_LIBRARIES(btkMathsLibraryTestsApp btkMathsLibrary ${CPPUNIT_LIBRARY} ${ITK_LIBRARIES})
    ADD_TEST(btkMathsLibraryTests btkMathsLibraryTestsApp)
ENDIF(BUILD_TESTING)

#---- Diffusion library -----------------------------------------------------------------

ADD_LIBRARY(btkDiffusionLibrary STATIC
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkMacro.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkFileHelper.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileHelper.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalDirection.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalHarmonics.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkGradientDirection.cxx ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkGradientDirection.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequence.cxx ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequence.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileReader.cxx ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileReader.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileWriter.cxx ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileWriter.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionTensorReconstructionFilter.cxx ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionTensorReconstructionFilter.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkSphericalHarmonicsDiffusionDecompositionFilter.cxx ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkSphericalHarmonicsDiffusionDecompositionFilter.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionModel.cxx ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionModel.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkTensorModel.cxx ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkTensorModel.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkOrientationDiffusionFunctionModel.cxx ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkOrientationDiffusionFunctionModel.h
)

TARGET_LINK_LIBRARIES(btkDiffusionLibrary btkMathsLibrary)
