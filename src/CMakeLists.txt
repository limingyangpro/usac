set(SRCEXT "cpp")
file(
    GLOB
    source_files
    *.${SRCEXT}
)


find_package(Eigen3 REQUIRED)
if(EIGEN3_INCLUDE_DIR)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_INCLUDE_DIR)

add_library( usaclib ${source_files})

