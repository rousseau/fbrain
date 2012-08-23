PROJECT(Tests)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Diffusion ${fbrain_SOURCE_DIR}/Code/InputOutput ${fbrain_SOURCE_DIR}/Code/Tools ${fbrain_SOURCE_DIR}/Code/Maths)

#---- VTK SandBox ----------------------------------------------------------------------------

ADD_EXECUTABLE(btkVTKSandBox btkVTKSandBox.cxx)
TARGET_LINK_LIBRARIES(btkVTKSandBox ${ITK_LIBRARIES} btkToolsLibrary btkMathsLibrary btkDiffusionLibrary vtkHybrid)

#---- Diffusion sequence ---------------------------------------------------------------------

ADD_EXECUTABLE(btkDiffusionSequenceTestApp ${fbrain_SOURCE_DIR}/Tests/btkDiffusionSequenceTest.cxx)
TARGET_LINK_LIBRARIES(btkDiffusionSequenceTestApp btkToolsLibrary btkMathsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})
ADD_TEST(btkDiffusionSequenceTest ${Tests_BINARY_DIR}/btkDiffusionSequenceTestApp ${Tests_SOURCE_DIR}/Resources/dwi.nii.gz ${Tests_BINARY_DIR}/dwi.nii.gz)

#---- Gradient directions --------------------------------------------------------------------

ADD_EXECUTABLE(btkGradientDirectionTestApp ${fbrain_SOURCE_DIR}/Tests/btkGradientDirectionTest.cxx)
TARGET_LINK_LIBRARIES(btkGradientDirectionTestApp btkToolsLibrary btkMathsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})
ADD_TEST(btkGradientDirectionTest ${Tests_BINARY_DIR}/btkGradientDirectionTestApp)

#---- Tractography ---------------------------------------------------------------------------

ADD_EXECUTABLE(btkTractographyTestApp btkTractographyTest.cxx)
TARGET_LINK_LIBRARIES(btkTractographyTestApp btkToolsLibrary btkMathsLibrary btkDiffusionLibrary ${ITK_LIBRARIES})
#ADD_TEST(btkTractographyTest ${Tests_BINARY_DIR/btkTractographyTestApp)
