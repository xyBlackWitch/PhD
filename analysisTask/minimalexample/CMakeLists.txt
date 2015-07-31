
set(INCLUDE_DIRECTORIES	
${ROOT_INCLUDE_DIR} 
${BASE_INCLUDE_DIRECTORIES}
${CMAKE_SOURCE_DIR}/rho
${CMAKE_SOURCE_DIR}/rho/RhoHistogram
${CMAKE_SOURCE_DIR}/rho/RhoBase
${CMAKE_SOURCE_DIR}/PndTools/AnalysisTools
${CMAKE_SOURCE_DIR}/AnalysisTask/AnalysisTask
)

include_directories(${INCLUDE_DIRECTORIES})
set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set(SRCS
AnalysisTask.cxx
)

set(LINKDEF anataskjennyLinkDef.h)
set(LIBRARY_NAME anataskjenny)

set(DEPENDENCIES Rho AnalysisTools)

GENERATE_LIBRARY()
