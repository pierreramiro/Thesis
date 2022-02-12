cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(thesis LANGUAGES CXX CUDA)

set_source_files_properties(src/main.cu PROPERTIES LANGUAGE "CUDA" )
add_executable(main src/main.cu)

set_target_properties(main PROPERTIES CUDA_SEPARABLE_COMPILATION ON)