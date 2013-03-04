PROJECT(Code)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#---- Directories -----------------------------------------------------------------------

SET(TOOLS_LIBRARY_SOURCE_DIR        ${Code_SOURCE_DIR}/Tools)
SET(MATHS_LIBRARY_SOURCE_DIR        ${Code_SOURCE_DIR}/Maths)
SET(DIFFUSION_LIBRARY_SOURCE_DIR    ${Code_SOURCE_DIR}/Diffusion)
SET(TRACTOGRAPHY_LIBRARY_SOURCE_DIR ${Code_SOURCE_DIR}/Tractography)

IF(BUILD_TESTING)
    SET(MATHS_TESTS_SOURCE_DIR ${MATHS_LIBRARY_SOURCE_DIR}/Tests)
ENDIF(BUILD_TESTING)

INCLUDE_DIRECTORIES(
    ${MATHS_LIBRARY_SOURCE_DIR}
    ${TOOLS_LIBRARY_SOURCE_DIR}
    ${DIFFUSION_LIBRARY_SOURCE_DIR}
    ${TRACTOGRAPHY_LIBRARY_SOURCE_DIR}
)

#---- Tools library ---------------------------------------------------------------------

SET(TOOLS_LIBRARY_HEADER
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkMacro.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkFileHelper.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkNrrdField.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkLandmarksFileReader.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkIOTransformHelper.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkImageHelper.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileHelper.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceHelper.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkCommandIterationUpdate.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkCommandProgressUpdate.h
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkIOImageHelper.h
)

SET(TOOLS_LIBRARY_SOURCES
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkFileHelper.cxx
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkNrrdField.cxx
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkLandmarksFileReader.cxx
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileHelper.cxx
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceHelper.cxx
    ${TOOLS_LIBRARY_SOURCE_DIR}/btkIOImageHelper.cxx
)

ADD_LIBRARY(btkToolsLibrary STATIC ${TOOLS_LIBRARY_HEADER} ${TOOLS_LIBRARY_SOURCES})
TARGET_LINK_LIBRARIES(btkToolsLibrary ${ITK_LIBRARIES})

#---- Maths library ---------------------------------------------------------------------

SET(MATHS_LIBRARY_HEADER
    ${MATHS_LIBRARY_SOURCE_DIR}/btkMathFunctions.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalDirection.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkLegendrePolynomial.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalHarmonics.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkMatrixOperations.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkHistogram.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkProbabilityDensity.h
    ${MATHS_LIBRARY_SOURCE_DIR}/btkNormalProbabilityDensity.h
)

SET(MATHS_LIBRARY_SOURCES
    ${MATHS_LIBRARY_SOURCE_DIR}/btkMathFunctions.cxx
    ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalDirection.cxx
    ${MATHS_LIBRARY_SOURCE_DIR}/btkLegendrePolynomial.cxx
    ${MATHS_LIBRARY_SOURCE_DIR}/btkSphericalHarmonics.cxx
    ${MATHS_LIBRARY_SOURCE_DIR}/btkMatrixOperations.cxx
    ${MATHS_LIBRARY_SOURCE_DIR}/btkHistogram.cxx
    ${MATHS_LIBRARY_SOURCE_DIR}/btkProbabilityDensity.cxx
    ${MATHS_LIBRARY_SOURCE_DIR}/btkNormalProbabilityDensity.cxx
)

ADD_LIBRARY(btkMathsLibrary STATIC ${MATHS_LIBRARY_HEADER} ${MATHS_LIBRARY_SOURCES})

IF(BUILD_TESTING)
    ADD_EXECUTABLE(btkMathsLibraryTestsApp ${MATHS_TESTS_SOURCE_DIR}/btkMathsLibraryTestsApp.cxx
        ${MATHS_TESTS_SOURCE_DIR}/btkMatrixOperationsTest.cxx
        ${MATHS_TESTS_SOURCE_DIR}/btkLegendrePolynomialTest.cxx
        ${MATHS_TESTS_SOURCE_DIR}/btkMathFunctionsTest.cxx
        ${MATHS_TESTS_SOURCE_DIR}/btkSphericalHarmonicsTest.cxx
        ${MATHS_TESTS_SOURCE_DIR}/btkNormalProbabilityDensityTest.cxx
    )
    TARGET_LINK_LIBRARIES(btkMathsLibraryTestsApp btkMathsLibrary ${CPPUNIT_LIBRARY} ${ITK_LIBRARIES})
    ADD_TEST(btkMathsLibraryTests btkMathsLibraryTestsApp)
ENDIF(BUILD_TESTING)

#---- Diffusion library -----------------------------------------------------------------

SET(DIFFUSION_LIBRARY_HEADER
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkGradientDirection.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequence.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileReader.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileWriter.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionGradientTable.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionTensorReconstructionFilter.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkSphericalHarmonicsDiffusionDecompositionFilter.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionModel.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkTensorModel.h
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkOrientationDiffusionFunctionModel.h
)

SET(DIFFUSION_LIBRARY_SOURCES
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkGradientDirection.cxx
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequence.cxx
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileReader.cxx
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionSequenceFileWriter.cxx
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionTensorReconstructionFilter.cxx
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkSphericalHarmonicsDiffusionDecompositionFilter.cxx
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkDiffusionModel.cxx
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkTensorModel.cxx
    ${DIFFUSION_LIBRARY_SOURCE_DIR}/btkOrientationDiffusionFunctionModel.cxx
)

ADD_LIBRARY(btkDiffusionLibrary STATIC ${DIFFUSION_LIBRARY_HEADER} ${DIFFUSION_LIBRARY_SOURCES})
TARGET_LINK_LIBRARIES(btkDiffusionLibrary btkToolsLibrary btkMathsLibrary)

#---- Tractography library --------------------------------------------------------------

SET(TRACTOGRAPHY_LIBRARY_HEADER
    ${TRACTOGRAPHY_LIBRARY_SOURCE_DIR}/btkTractographyAlgorithm.h
    ${TRACTOGRAPHY_LIBRARY_SOURCE_DIR}/btkStreamlineTractographyAlgorithm.h
    ${TRACTOGRAPHY_LIBRARY_SOURCE_DIR}/btkParticleFilteringTractographyAlgorithm.h
    ${TRACTOGRAPHY_LIBRARY_SOURCE_DIR}/btkPolyDataColorLinesByOrientation.h
)

SET(TRACTOGRAPHY_LIBRARY_SOURCES
    ${TRACTOGRAPHY_LIBRARY_SOURCE_DIR}/btkTractographyAlgorithm.cxx
    ${TRACTOGRAPHY_LIBRARY_SOURCE_DIR}/btkStreamlineTractographyAlgorithm.cxx
    ${TRACTOGRAPHY_LIBRARY_SOURCE_DIR}/btkParticleFilteringTractographyAlgorithm.cxx
    ${TRACTOGRAPHY_LIBRARY_SOURCE_DIR}/btkPolyDataColorLinesByOrientation.cxx
)

ADD_LIBRARY(btkTractographyLibrary STATIC ${TRACTOGRAPHY_LIBRARY_HEADER} ${TRACTOGRAPHY_LIBRARY_SOURCES})
TARGET_LINK_LIBRARIES(btkTractographyLibrary ${ITK_LIBRARIES})
