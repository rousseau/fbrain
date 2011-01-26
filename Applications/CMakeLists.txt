PROJECT(Applications)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INCLUDE_DIRECTORIES(${fbrain_SOURCE_DIR}/Code/Reconstruction)

ADD_EXECUTABLE(btkNLMDenoising btkNLMDenoising.cxx )
TARGET_LINK_LIBRARIES(btkNLMDenoising ITKIO)

ADD_EXECUTABLE(btkImageReconstruction btkImageReconstruction.cxx )
TARGET_LINK_LIBRARIES(btkImageReconstruction ITKIO ITKNumerics)


