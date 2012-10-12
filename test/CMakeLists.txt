# Libraries used by all builds
SET(common_link_libraries ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${FFTW_LIBRARIES})

# Loop for each test and create an executable and create a test that
# runs the executable.
#SET(test_names peak_search sss_detect tfg xcorr_pss)
#FOREACH (TN ${test_names})
#  ADD_EXECUTABLE(test_${TN} test_${TN}.cpp)
#  TARGET_LINK_LIBRARIES (test_${TN} general LTE_MISC curses)
#  TARGET_LINK_LIBRARIES (test_${TN} debug itpp_debug ${common_link_libraries})
#  TARGET_LINK_LIBRARIES (test_${TN} optimized itpp ${common_link_libraries})
#  ADD_TEST (${TN} test_${TN} "${PROJECT_SOURCE_DIR}/test/test_${TN}.it")
#  SET_TESTS_PROPERTIES(${TN} PROPERTIES PASS_REGULAR_EXPRESSION passed)
#ENDFOREACH (TN)

