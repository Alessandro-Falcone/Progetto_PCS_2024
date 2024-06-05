list(APPEND srcUt_headers "${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.hpp")
list(APPEND srcUt_sources "${CMAKE_CURRENT_SOURCE_DIR}/Utils2poligoni.cpp")
list(APPEND srcUt_sources "${CMAKE_CURRENT_SOURCE_DIR}/Utils2poligoni.hpp")


set(srcUt_headers ${srcUt_headers} PARENT_SCOPE)
set(srcUt_sources ${srcUt_sources} PARENT_SCOPE)
set(srcUt_headers ${srcUt_headers} PARENT_SCOPE)
set(srcUt_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
