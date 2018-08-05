# Adding all the test suites
add_executable(test_baseline test_JPEGDecoderBaseline.cpp)
target_link_libraries(test_baseline gtest_main)
add_test(NAME all_test_baseline COMMAND test_baseline)

add_executable(test_huffman test_JPEGDecoderHuffman.cpp)
target_link_libraries(test_huffman gtest_main)
add_test(NAME all_test_huffman COMMAND test_huffman)

add_executable(test_marker test_JPEGDecoderMarker.cpp)
target_link_libraries(test_marker gtest_main)
add_test(NAME all_test_marker COMMAND test_marker)

add_executable(test_parsing test_JPEGDecoderParsing.cpp)
target_link_libraries(test_parsing gtest_main)
add_test(NAME all_test_parsing COMMAND test_parsing)

add_executable(test_utility test_JPEGDecoderUtility.cpp)
target_link_libraries(test_utility gtest_main)
add_test(NAME all_test_utility COMMAND test_utility)