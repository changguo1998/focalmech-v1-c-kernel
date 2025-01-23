cmake_minimum_required(VERSION 3.22)
project(focalmech_C LANGUAGES C)

include_directories(include)

find_package(OpenMP)

add_executable(fm include/fm.h include/kernel.h include/types.h src/fm.c src/kernel.c src/types.c
        src/io.c include/io.h)
target_include_directories(fm PRIVATE include)
target_link_libraries(fm PRIVATE m)
install(TARGETS fm)

add_executable(fm_omp include/fm.h include/kernel.h include/types.h src/fm_omp.c src/kernel.c src/types.c
        src/io.c include/io.h)
target_include_directories(fm_omp PRIVATE include)
target_link_libraries(fm_omp PRIVATE m OpenMP::OpenMP_C)
install(TARGETS fm_omp)

#add_executable(arrayio test/array_io.c)
#install(TARGETS arrayio)
