cmake_minimum_required(VERSION 3.10)
project(SolveQP C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

include_directories(include)

set(SOURCE_FILES
        src/solve_qp.cpp
        src/solve.QP.c
        src/qr-update.c
        src/linear-algebra.c
        src/decomposeQP.cpp
        src/findSigExposures.cpp
        main.cpp
)

add_executable(solve_qp ${SOURCE_FILES})