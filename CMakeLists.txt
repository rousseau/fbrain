PROJECT(fbrain)

cmake_minimum_required(VERSION 2.6)
OPTION(USE_OMP "Use OpenMP" ON)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp  ")

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

SUBDIRS( 
  Applications
  Utilities
)


