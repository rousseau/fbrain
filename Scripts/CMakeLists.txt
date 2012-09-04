PROJECT(Scripts)
INCLUDE_REGULAR_EXPRESSION("^.*$")

INSTALL(FILES
    ${Scripts_SOURCE_DIR}/btkPrepareAtlasData.py
    ${Scripts_SOURCE_DIR}/btkCreateTemplate.py
DESTINATION bin)
