cmake_minimum_required(VERSION 3.5)
project(gpsp)

include_directories("/Library/cplex1263/cplex/include")
include_directories("/Library/cplex1263/cplex/include/ilcplex")

set(CMAKE_C_COMPILER "/usr/local/bin/clang-omp")
set(CMAKE_CXX_COMPILER "/usr/local/bin/clang-omp++")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto -fopenmp -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto -lm -fopenmp -DCPX")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/Library/cplex1263/cplex/lib/x86-64_osx/static_pic -L/Library/cplex1263/concert/lib/x86-64_osx/static_pic")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto -lm -fopenmp -lpthread -lcplex")

set(SOURCE_FILES
        dict_int.c
        dict_int.h
        instance.c
        instance.h
        lahc.c
        lahc.h
        learning.c
        learning.h
        list_int.c
        list_int.h
        long_compl_path.c
        long_compl_path.h
        lp.cpp
        lp.h
        lslahc.c
        macros.h
        memory.c
        memory.h
        mode_set.c
        mode_set.h
        ms_solver_mip.c
        ms_solver_mip.h
        neighborhood.c
        neighborhood.h
        node_heap.c
        node_heap.h
        rrusage.c
        rrusage.h
        solution.c
        solution.h
        str_utils.c
        str_utils.h
        test.c
        test.h
        tokenizer.c
        tokenizer.h
        top_sort.c
        top_sort.h
        vec_char.c
        vec_char.h
        vec_double.c
        vec_double.h
        vec_int.c
        vec_int.h
        vec_str.c
        vec_str.h)

add_executable(src ${SOURCE_FILES})