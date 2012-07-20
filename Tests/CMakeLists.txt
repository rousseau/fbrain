PROJECT(Tests)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Diffusion ${fbrain_SOURCE_DIR}/Code/InputOutput ${fbrain_SOURCE_DIR}/Code/Tools ${fbrain_SOURCE_DIR}/Code/Maths)

#---- Diffusion sequence ---------------------------------------------------------------------

ADD_EXECUTABLE(btkDiffusionSequenceTestApp ${fbrain_SOURCE_DIR}/Tests/btkDiffusionSequenceTest.cxx
    ${fbrain_SOURCE_DIR}/Code/Tools/btkFileHelper.cxx ${fbrain_SOURCE_DIR}/Code/Tools/btkFileHelper.h
    ${fbrain_SOURCE_DIR}/Code/Tools/btkDiffusionSequenceHelper.h
    ${fbrain_SOURCE_DIR}/Code/Tools/btkDiffusionSequenceFileHelper.cxx ${fbrain_SOURCE_DIR}/Code/Tools/btkDiffusionSequenceFileHelper.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkGradientDirection.cxx ${fbrain_SOURCE_DIR}/Code/Diffusion/btkGradientDirection.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequence.cxx ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequence.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequenceFileReader.cxx ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequenceFileReader.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequenceFileWriter.cxx ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequenceFileWriter.h
)
TARGET_LINK_LIBRARIES(btkDiffusionSequenceTestApp ${ITK_LIBRARIES})
ADD_TEST(btkDiffusionSequenceTest ${Tests_BINARY_DIR}/btkDiffusionSequenceTestApp ${Tests_SOURCE_DIR}/Resources/dwi.nii.gz ${Tests_BINARY_DIR}/dwi.nii.gz)

#---- Gradient directions --------------------------------------------------------------------

ADD_EXECUTABLE(btkGradientDirectionTestApp ${fbrain_SOURCE_DIR}/Tests/btkGradientDirectionTest.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkGradientDirection.cxx ${fbrain_SOURCE_DIR}/Code/Diffusion/btkGradientDirection.h
)
TARGET_LINK_LIBRARIES(btkGradientDirectionTestApp ${ITK_LIBRARIES})
ADD_TEST(btkGradientDirectionTest ${Tests_BINARY_DIR}/btkGradientDirectionTestApp)

#---- Tractography ---------------------------------------------------------------------------

ADD_EXECUTABLE(btkTractographyTestApp btkTractographyTest.cxx
    ${fbrain_SOURCE_DIR}/Code/Maths/btkSphericalHarmonics.h ${fbrain_SOURCE_DIR}/Code/Maths/btkSphericalHarmonics.cxx
    ${fbrain_SOURCE_DIR}/Code/Maths/btkLegendrePolynomial.h ${fbrain_SOURCE_DIR}/Code/Maths/btkLegendrePolynomial.cxx
    ${fbrain_SOURCE_DIR}/Code/Maths/btkMathFunctions.h ${fbrain_SOURCE_DIR}/Code/Maths/btkMathFunctions.cxx
    ${fbrain_SOURCE_DIR}/Code/Maths/btkSphericalDirection.h ${fbrain_SOURCE_DIR}/Code/Maths/btkSphericalDirection.cxx
    ${fbrain_SOURCE_DIR}/Code/Tools/btkFileHelper.cxx ${fbrain_SOURCE_DIR}/Code/Tools/btkFileHelper.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionModel.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkSphericalHarmonicsDiffusionDecomposition.h ${fbrain_SOURCE_DIR}/Code/Diffusion/btkSphericalHarmonicsDiffusionDecomposition.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionTensor.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkGradientDirection.h ${fbrain_SOURCE_DIR}/Code/Diffusion/btkGradientDirection.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequence.h ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequence.cxx
    ${fbrain_SOURCE_DIR}/Code/Tools/btkDiffusionSequenceFileHelper.cxx ${fbrain_SOURCE_DIR}/Code/Tools/btkDiffusionSequenceFileHelper.h
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequenceFileReader.h ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequenceFileReader.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequenceFileWriter.h ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionSequenceFileWriter.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionTensorReconstructionFilter.h ${fbrain_SOURCE_DIR}/Code/Diffusion/btkDiffusionTensorReconstructionFilter.cxx
    ${fbrain_SOURCE_DIR}/Code/Diffusion/btkSphericalHarmonicsDiffusionDecompositionReconstructionFilter.h ${fbrain_SOURCE_DIR}/Code/Diffusion/btkSphericalHarmonicsDiffusionDecompositionReconstructionFilter.cxx
)
TARGET_LINK_LIBRARIES(btkTractographyTestApp ${ITK_LIBRARIES})
#ADD_TEST(btkTractographyTest ${Tests_BINARY_DIR/btkTractographyTestApp)
