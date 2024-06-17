list(APPEND srcParte2_headers "${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.hpp")
list(APPEND srcParte2_sources "${CMAKE_CURRENT_SOURCE_DIR}/Utils2poligoni.cpp")
list(APPEND srcParte2_headers "${CMAKE_CURRENT_SOURCE_DIR}/Utils2poligoni.hpp")


set(srcParte2_headers ${srcParte2_headers} PARENT_SCOPE)
set(srcParte2_sources ${srcParte2_sources} PARENT_SCOPE)
set(srcParte2_headers ${srcParte2_headers} PARENT_SCOPE)
set(srcParte2_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
