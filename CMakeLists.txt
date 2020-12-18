cmake_minimum_required(VERSION 3.16)
project(prime-gap)

set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_LIST_DIR}/cmake
    ${CMAKE_MODULE_PATH}
)

find_package(GMP REQUIRED)
message(STATUS "GMP_INCLUDES: ${GMP_INCLUDES} ${GMP_LIBRARIES}")
link_libraries(${GMP_LIBRARIES})

# When gmp 6.3 is out with mpz_prevprime support and in ubuntu
#link_libraries(gmp)

set(CMAKE_CXX_STANDARD 14)

add_executable(combined_sieve   combined_sieve.cpp gap_common.cpp modulo_search.cpp)
add_executable(gap_stats        gap_stats.cpp gap_common.cpp)
add_executable(gap_test_simple  gap_test_simple.cpp gap_common.cpp)
add_executable(benchmark        misc/benchmark.cpp modulo_search.cpp)

target_link_libraries(combined_sieve sqlite3)
target_link_libraries(gap_stats sqlite3)
target_link_libraries(gap_test_simple sqlite3)


target_link_libraries(combined_sieve primesieve)
target_link_libraries(benchmark primesieve)
