cmake_minimum_required(VERSION 3.16)
project(NeoEuclidSearch VERSION 0.1.0 LANGUAGES CXX)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(neo_euclid
    src/engine/sample_batch.cpp
    src/geometry.cpp
    src/state.cpp
    src/search.cpp
    src/problems/fixed_angle.cpp
    src/problems/fixed_angle_mitm.cpp
    src/problems/obtuse.cpp
    src/problems/parabola.cpp
    src/problems/circumcircle.cpp
)
target_include_directories(neo_euclid PUBLIC include)
target_link_libraries(neo_euclid PUBLIC Threads::Threads)
if(MSVC)
    target_compile_options(neo_euclid PRIVATE /W4 /permissive-)
else()
    target_compile_options(neo_euclid PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(neo-euclid src/main.cpp)
target_link_libraries(neo-euclid PRIVATE neo_euclid)

include(CTest)
if(BUILD_TESTING)
    add_executable(neo_euclid_tests tests/test_main.cpp)
    target_link_libraries(neo_euclid_tests PRIVATE neo_euclid)
    add_test(NAME neo_euclid_tests COMMAND neo_euclid_tests)
endif()
