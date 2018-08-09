# Adding all the test suites

add_executable(test_huffman test_JPEGHuffman.cpp)
target_link_libraries(test_huffman gmock_main JPEGDecoder  ${OpenCV_LIBS})
add_test(NAME all_test_huffman COMMAND test_huffman)

add_executable(test_parsing test_JPEGParsing.cpp)
target_link_libraries(test_parsing gmock_main JPEGDecoder  ${OpenCV_LIBS})
add_test(NAME all_test_parsing COMMAND test_parsing)

add_executable(test_utility test_JPEGUtility.cpp )
target_link_libraries(test_utility gmock_main JPEGDecoder  ${OpenCV_LIBS})
add_test(NAME all_test_utility COMMAND test_utility)
