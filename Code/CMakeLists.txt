PROJECT(Code)
INCLUDE_REGULAR_EXPRESSION("^.*$")

#---- Directories -----------------------------------------------------------------------

SET(TOOLS_LIBRARY_SOURCE_DIR ${Code_SOURCE_DIR}/Tools)
SET(MATHS_LIBRARY_SOURCE_DIR ${Code_SOURCE_DIR}/Maths)

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
#TARGET_LINK_LIBRARIES
