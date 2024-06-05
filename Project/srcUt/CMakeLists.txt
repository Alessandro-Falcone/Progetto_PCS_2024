list(APPEND srcUt_headers "${CMAKE_CURRENT_SOURCE_DIR}/DFN.hpp")
list(APPEND srcUt_sources "${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp")
list(APPEND srcUt_headers "${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp")


set(srcUt_headers ${srcUt_headers} PARENT_SCOPE)
set(srcUt_sources ${srcUt_sources} PARENT_SCOPE)
set(srcUt_headers ${srcUt_headers} PARENT_SCOPE)
set(srcUt_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
