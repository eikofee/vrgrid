cmake_minimum_required(VERSION 3.10)

add_executable(vrgrid vrgrid.cpp)
INCLUDE_DIRECTORIES(
    ${WulipBaseInclude}
    ${WulipCoreInclude}
    ${WulipAlgoInclude}
    ${THIRDPARTY_DIR}/freetype-gl
    ${THIRDPARTY_DIR}/optional
    ${WulipRendererInclude}

)

FIND_PACKAGE(GLUT REQUIRED)
TARGET_LINK_LIBRARIES(vrgrid
    pthread
    ${LibWulipAlgoName}
    ${LibWulipRendererName}
    ${GLUT_LIBRARY}
)

INSTALL(TARGETS vrgrid  RUNTIME DESTINATION ${WulipBinInstallDir} COMPONENT demo)
#resources installation
  INSTALL(
    DIRECTORY ${PROJECT_SOURCE_DIR}/renderer/textures ${PROJECT_SOURCE_DIR}/renderer/shaders ${PROJECT_SOURCE_DIR}/share/fonts DESTINATION ${WulipBinInstallDir} COMPONENT demo)
