PROJECT(Viewer)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INCLUDE_DIRECTORIES(
    ${TCLAP_DIRECTORY}
    ${Code_SOURCE_DIR}/Diffusion
    ${Code_SOURCE_DIR}/ImageFilters
    ${Code_SOURCE_DIR}/Tools
    ${Code_SOURCE_DIR}/Maths
    ${Code_SOURCE_DIR}/Transformations
    ${Code_SOURCE_DIR}/Visualization
)

ADD_EXECUTABLE( btkSlicesMotionViewer btkSlicesMotionViewer.cxx
${Code_SOURCE_DIR}/Visualization/btkRenderPolydata.cxx
${Code_SOURCE_DIR}/Visualization/btkRenderPolydata.h)
TARGET_LINK_LIBRARIES(btkSlicesMotionViewer btkMathsLibrary ${ITK_LIBRARIES} vtkHybrid)

ADD_EXECUTABLE(btkViewer btkViewer.cxx)
TARGET_LINK_LIBRARIES(btkViewer btkToolsLibrary btkDiffusionLibrary ${ITK_LIBARIES} vtkHybrid vtkWidgets)

ADD_EXECUTABLE(btkMRIViewer btkMRIViewer.cxx)
TARGET_LINK_LIBRARIES(btkMRIViewer btkToolsLibrary ${ITK_LIBARIES} vtkHybrid vtkWidgets)

INSTALL(TARGETS btkSlicesMotionViewer
        btkViewer
        btkMRIViewer
        DESTINATION bin)

