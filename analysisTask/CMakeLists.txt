
set(INCLUDE_DIRECTORIES	
${ROOT_INCLUDE_DIR} 
${BASE_INCLUDE_DIRECTORIES}
${CMAKE_SOURCE_DIR}/rho
${CMAKE_SOURCE_DIR}/rho/RhoHistogram
${CMAKE_SOURCE_DIR}/rho/RhoBase
${CMAKE_SOURCE_DIR}/rho/RhoSelector
${CMAKE_SOURCE_DIR}/PndTools/AnalysisTools
${CMAKE_SOURCE_DIR}/PndTools/AnalysisTools/Fitter
${CMAKE_SOURCE_DIR}/AnalysisTask
${CMAKE_SOURCE_DIR}/pnddata
${CMAKE_SOURCE_DIR}/pnddata/PidData
)

include_directories(${INCLUDE_DIRECTORIES})
set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set(SRCS
AnalysisTask.cxx
AnalysisTaskLambda0.cxx
AnalysisTaskXi1820.cxx
AnalysisTaskAntiXi1820.cxx
)

set(LINKDEF anataskjennyLinkDef.h)
set(LIBRARY_NAME anataskjenny)

set(DEPENDENCIES PndData Rho AnalysisTools)

GENERATE_LIBRARY()
