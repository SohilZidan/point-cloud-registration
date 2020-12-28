cmake_minimum_required(VERSION 3.11)
project( ICP )
# Find nanoflannConfig.cmake:
find_package(nanoflann)
find_package(Eigen3)


# set the c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(testTrICP src/main.cpp happly.h src/gsec.cpp src/utils.cpp src/icp.cpp)

target_include_directories(testTrICP PRIVATE ${EIGEN3_INCLUDE_DIR})

# Make sure the include path is used:
target_link_libraries(testTrICP nanoflann::nanoflann)
