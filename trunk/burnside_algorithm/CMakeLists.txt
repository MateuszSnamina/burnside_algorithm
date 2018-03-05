cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
file(GLOB sources src/*.cpp)
#message("${sources}")
add_executable(burnside ${sources})
target_include_directories(burnside PUBLIC include)
target_link_libraries(burnside armadillo)
