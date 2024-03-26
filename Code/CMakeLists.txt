cmake_minimum_required(VERSION 3.13)
project(lokimc)

set(CMAKE_CXX_STANDARD 17)

# --- Compiler flags --- #

# openmp package
FIND_PACKAGE(OpenMP REQUIRED)

if(MSVC)
    message(STATUS "Using MSVC compiler")

    # #define _USE_MATH_DEFINES for M_PI and M_E
    # #define _CRT_SECURE_NO_DEPRECATE to surpress fopen warnings (fopen_s is not available on other platforms)
    # #define _MSVC to alter a regular expression in the Parse class (since MSVC treats the $ command differently)
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_DEPRECATE -D_MSVC)

    # warning flag: /W4
    # there is no alternative for -march=native in MSVC, therefore important relevant compiler flags have to be set manually for now (avx, avx2, fma etc.).
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG /openmp ${OpenMP_C_FLAGS}")
 
    # set necessary variables in MSVC
    set(Boost_USE_STATIC_LIBS OFF) 
    set(Boost_USE_MULTITHREADED ON)  
    set(Boost_USE_STATIC_RUNTIME OFF) 

else() # clang or g++
    # warning flags: -Wall -Wextra -Werror
    message(STATUS "Using g++ or Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 ${OpenMP_C_FLAGS}") 
    # turning off deprecation warnings
    set(DEPRECATION_FLAG "-Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${DEPRECATION_FLAG}")    
endif()

message(${CMAKE_CXX_FLAGS_RELEASE})

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# --- loki-mc --- #
include_directories(./)
AUX_SOURCE_DIRECTORY(LoKI-MC/Sources srcFiles1)
AUX_SOURCE_DIRECTORY(External/MathParser srcFiles2)
add_executable(lokimc ${srcFiles1} ${srcFiles2})

# Boost library
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# GSL library
find_package(GSL REQUIRED) 
target_link_libraries(lokimc GSL::gsl GSL::gslcblas)

# OpenMP library
find_package(OpenMP REQUIRED)
target_link_libraries(lokimc OpenMP::OpenMP_CXX)