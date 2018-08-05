cmake_minimum_required(VERSION 3.1.0)
project(jpeg_decoder VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)
set(Tutorial_VERSION_MAJOR 0)
set(Tutorial_VERSION_MINOR 1)

include(CTest)
enable_testing()

include_directories("include/")
find_package( OpenCV REQUIRED )



#-------------------------------------------------
#                GOOGLE TEST CONF
#-------------------------------------------------

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )

if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                 ${CMAKE_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

# Adding all the test suites
add_executable(test_baseline test/test_JPEGDecoderBaseline.cpp)
target_link_libraries(test_baseline gtest_main)
add_test(NAME all_test_baseline COMMAND test_baseline)

add_executable(test_huffman test/test_JPEGDecoderHuffman.cpp)
target_link_libraries(test_huffman gtest_main)
add_test(NAME all_test_huffman COMMAND test_huffman)

add_executable(test_marker test/test_JPEGDecoderMarker.cpp)
target_link_libraries(test_marker gtest_main)
add_test(NAME all_test_marker COMMAND test_marker)

add_executable(test_parsing test/test_JPEGDecoderParsing.cpp)
target_link_libraries(test_parsing gtest_main)
add_test(NAME all_test_parsing COMMAND test_parsing)

add_executable(test_utility test/test_JPEGDecoderUtility.cpp)
target_link_libraries(test_utility gtest_main)
add_test(NAME all_test_utility COMMAND test_utility)

# Adding the main code executable
add_executable(jpeg_decoder src/main.cpp src/JPEGDecoder.cpp include/JPEGDecoder.hpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_link_libraries(jpeg_decoder ${OpenCV_LIBS})