PROJECT(Scripts)
INCLUDE_REGULAR_EXPRESSION("^.*$")

IF(ENABLE_SCRIPT)
    INSTALL(FILES
        ${Scripts_SOURCE_DIR}/btkAffineTools.py
        ${Scripts_SOURCE_DIR}/btkPatientsTools.py
        ${Scripts_SOURCE_DIR}/btkGaussianKernel.py
        ${Scripts_SOURCE_DIR}/btkPrepareAtlasData.py
        ${Scripts_SOURCE_DIR}/btkCreateTemplate.py
        ${Scripts_SOURCE_DIR}/btkCreateLongitudinalAtlas.py
        ${Scripts_SOURCE_DIR}/btkMaskExtractionProcess.py
        ${Scripts_SOURCE_DIR}/btkMaskExtraction.py
    DESTINATION bin
    PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
ENDIF(ENABLE_SCRIPT)
