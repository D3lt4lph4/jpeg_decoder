cmake_minimum_required(VERSION 3.1.0)
project(jpeg_decoder VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)
set(Tutorial_VERSION_MAJOR 0)
set(Tutorial_VERSION_MINOR 1)
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG")
include(CTest)
enable_testing()

if(WIN32)
  set(OpenCV_DIR "C:/Program Files (x86)/opencv/build")
endif()

find_package(OpenCV REQUIRED)
find_package(Boost COMPONENTS thread system filesystem log REQUIRED) 

message(STATUS "Found opencv library : ${OpenCV_LIBS}\n")
message(STATUS "Found boost library : ${Boost_LIBRARIES}\n")


#-------------------------------------------------
#                SWIG CONF
#-------------------------------------------------

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs 2.7 REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(decoder.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(decoder.i PROPERTIES SWIG_FLAGS "-includeall")
file(GLOB_RECURSE Swig_SRC "src/*.cpp")

SWIG_ADD_MODULE(decoder python decoder.i ${Swig_SRC})
SWIG_LINK_LIBRARIES(decoder ${OpenCV_LIBS} ${Boost_LIBRARIES})
SET_SOURCE_FILES_PROPERTIES(${swig_generated_file_fullname} PROPERTIES COMPILE_FLAGS "-fPIC")
SWIG_LINK_LIBRARIES(decoder ${PYTHON_LIBRARIES})

#-------------------------------------------------
#               END SWIG CONF
#-------------------------------------------------

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
  include_directories(BEFORE SYSTEM
    "${gtest_SOURCE_DIR}/include" "${gmock_SOURCE_DIR}/include")
else()
  target_include_directories(gmock_main SYSTEM BEFORE INTERFACE
    "${gtest_SOURCE_DIR}/include" "${gmock_SOURCE_DIR}/include")
endif()

#-------------------------------------------------
#               END GOOGLE TEST CONF
#-------------------------------------------------

include_directories("include/")

add_subdirectory(test)

# Linking boost library (dynamic)
add_definitions(-DBOOST_LOG_DYN_LINK)

# Creating library for the decoder
add_library(JPEGDecoder SHARED src/JPEGDecoder.cpp src/JPEGHuffmanDecoder.cpp src/JPEGUtility.cpp src/JPEGParser.cpp)
target_link_libraries(JPEGDecoder ${Boost_LIBRARIES})

# Adding the main code executable
add_executable(jpeg_decoder src/main.cpp)

target_link_libraries(jpeg_decoder ${OpenCV_LIBS} ${Boost_LIBRARIES} JPEGDecoder)