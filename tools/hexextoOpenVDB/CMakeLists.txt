cmake_minimum_required(VERSION 3.1)
project(hexextoOpenVDB)

find_package(OpenMP REQUIRED)

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)	  

add_executable(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} ../../src)

# Use C++14
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} cubecover polyscope OpenVDB::openvdb OpenMP::OpenMP_CXX)
