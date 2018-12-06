cmake_minimum_required(VERSION 3.1.0)
project(jpeg_decoder VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-g")
include(CTest)
enable_testing()

# We need boost for the logs (might be removed)
find_package(Boost COMPONENTS thread system filesystem log REQUIRED)
message(STATUS "Found boost library : ${Boost_LIBRARIES}\n")
add_definitions(-DBOOST_LOG_DYN_LINK)

# All the hpp files
include_directories("include/")

# Creating library for the decoder
add_library(JPEGDecoder SHARED src/JPEGDecoder.cpp src/JPEGHuffmanDecoder.cpp src/JPEGUtility.cpp src/JPEGParser.cpp src/JPEGUtils.cpp)
target_link_libraries(JPEGDecoder ${Boost_LIBRARIES})

# If we are not in debug mode then compile for deployement
IF( (NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES RELEASE))
  # Part for the python wrapper
  add_subdirectory(pybind11)

  # pybind11 package creation
  pybind11_add_module(jpegdecoder src/main_python.cpp src/JPEGUtils_python.cpp src/JPEGUtility_python.cpp src/JPEGDecoder_python.cpp)
  target_link_libraries(jpegdecoder PRIVATE ${Boost_LIBRARIES} JPEGDecoder)

ELSEIF( CMAKE_BUILD_TYPE MATCHES DEBUG )
  set(CMAKE_CXX_FLAGS "-g")
  # We set opencv only if we are in debug
  IF( WIN32 )
    set(OpenCV_DIR "C:/Program Files (x86)/opencv/build")
  ENDIF( WIN32 )
  set(CMAKE_BUILD_TYPE Debug)
  find_package(OpenCV REQUIRED)
  message(STATUS "Found opencv library : ${OpenCV_LIBS}\n")
  
  # Setting the flag for the DEBUG part of the code
  set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG")

  # Initializing the test files for the debug
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

  # Adding the test subdirectory
  add_subdirectory(test)

  # Finally we create a main program that can be used to use the library.
  add_executable(jpeg_decoder src/main.cpp)

  target_link_libraries(jpeg_decoder ${OpenCV_LIBS} ${Boost_LIBRARIES} JPEGDecoder)
ELSE( (NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES RELEASE) )
  message("Error, the build should be either in debug or release mode.")
ENDIF( (NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES RELEASE) )
