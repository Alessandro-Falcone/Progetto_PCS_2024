list(APPEND srcGoogleTest_headers "${CMAKE_CURRENT_SOURCE_DIR}/GoogleTest.hpp")


set(srcGoogleTest_sources ${srcGoogleTest_sources} PARENT_SCOPE)
set(srcGoogleTest_headers ${srcGoogleTest_headers} PARENT_SCOPE)
set(srcGoogleTest_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
