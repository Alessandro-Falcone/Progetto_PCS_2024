list(APPEND srcParaview_sources "${CMAKE_CURRENT_SOURCE_DIR}/Paraview.cpp")
list(APPEND srcParaview_headers "${CMAKE_CURRENT_SOURCE_DIR}/Paraview.hpp")


set(srcParaview_sources ${srcParaview_sources} PARENT_SCOPE)
set(srcParaview_headers ${srcParaview_headers} PARENT_SCOPE)
set(srcParaview_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
