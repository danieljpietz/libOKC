cmake_minimum_required(VERSION 3.23)
project(CPP_NE_2)

set(CMAKE_CXX_STANDARD 20)

include_directories(lib src)

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-march=native" _march_native_works)
check_cxx_compiler_flag("-xHost" _xhost_works)

set(_CXX_FLAGS)
if(_march_native_works)
    message(STATUS "Using processor's vector instructions (-march=native compiler flag set)")
    set(_CXX_FLAGS "-march=native")
elseif(_xhost_works)
    message(STATUS "Using processor's vector instructions (-xHost compiler flag set)")
    set(_CXX_FLAGS "-xHost")
else()
    message(STATUS "No suitable compiler flag found for vectorization")
endif()

add_executable(CPP_NE_2 main.cpp src/okc_types.h src/okc.h src/okc_system.h src/okc_math.h src/okc_link.h src/okc_differential_types.h)


#target_compile_options(CPP_NE_2 PRIVATE ${_CXX_FLAGS})
