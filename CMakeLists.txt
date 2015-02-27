cmake_minimum_required(VERSION 2.8)
project(SIMBODY_ROLLING_DISC)

find_package(Simbody REQUIRED)

include_directories(${Simbody_INCLUDE_DIR})
link_directories(${Simbody_LIB_DIR})

set(CMAKE_CXX_FLAGS "-Wall -std=c++11 ${CMAKE_CXX_FLAGS}")

add_executable(rolling_disc src/disc.cc)
target_link_libraries(rolling_disc ${Simbody_LIBRARIES})

