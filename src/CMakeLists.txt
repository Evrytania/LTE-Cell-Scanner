# Create a library of all the shared functions.
add_library(LTE_MISC capbuf.cpp constants.cpp itpp_ext.cpp macros.cpp searcher.cpp common.cpp dsp.cpp lte_lib.cpp from_osmocom.cpp)

SET (common_link_libs ${Boost_LIBRARIES} ${Boost_THREAD_LIBRARY} ${LAPACK_LIBRARIES} ${FFTW_LIBRARIES} ${RTLSDR_LIBRARIES} ${CURSES_LIBRARIES})

# Create the cell search executable
ADD_EXECUTABLE (CellSearch CellSearch.cpp)
TARGET_LINK_LIBRARIES (CellSearch general LTE_MISC)
TARGET_LINK_LIBRARIES (CellSearch debug itpp_debug ${common_link_libs})
TARGET_LINK_LIBRARIES (CellSearch optimized itpp ${common_link_libs})

# Create the cell tracker
ADD_EXECUTABLE (LTE-Tracker LTE-Tracker.cpp producer_thread.cpp tracker_thread.cpp searcher_thread.cpp display_thread.cpp)
TARGET_LINK_LIBRARIES (LTE-Tracker general LTE_MISC)
TARGET_LINK_LIBRARIES (LTE-Tracker debug itpp_debug ${common_link_libs})
TARGET_LINK_LIBRARIES (LTE-Tracker optimized itpp ${common_link_libs})

# Code use to test whether the rtl_sdr dongle is dropping samples
# or not.
#ADD_EXECUTABLE (rtl_sdr_check rtl_sdr_check.cpp)
#TARGET_LINK_LIBRARIES (rtl_sdr_check general LTE_MISC curses)
#TARGET_LINK_LIBRARIES (rtl_sdr_check debug itpp_debug ${common_link_libs})
#TARGET_LINK_LIBRARIES (rtl_sdr_check optimized itpp ${common_link_libs})

# Experimental code...
#ADD_EXECUTABLE (exp exp.cpp)
#TARGET_LINK_LIBRARIES (exp general LTE_MISC)
#TARGET_LINK_LIBRARIES (exp debug itpp_debug ${common_link_libs})
#TARGET_LINK_LIBRARIES (exp optimized itpp ${common_link_libs})

# Create a test that tests the overall functionality of the cell searcher.
# The tests in the 'tests' subdirectory test some of the smaller internal
# functions.
#ADD_TEST(FullTest CellSearch -s 739000000 -l -d "${PROJECT_SOURCE_DIR}/test")
#SET_TESTS_PROPERTIES(FullTest PROPERTIES PASS_REGULAR_EXPRESSION cell.ID..271)

INSTALL( TARGETS CellSearch LTE-Tracker DESTINATION bin )

