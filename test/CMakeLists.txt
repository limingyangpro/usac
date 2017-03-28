find_package ( OpenCV REQUIRED )
find_package( Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
find_package( Eigen3 REQUIRED)

set(SRCEXT "cpp")
file(
    GLOB
    source_files
    *.${SRCEXT}
)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 

add_executable( runAllTests ${source_files})

enable_testing()
ADD_TEST (runAllTests runAllTests)

target_link_libraries( runAllTests 
                    usaclib 
                    ${OpenCV_LIBS})

