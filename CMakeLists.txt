cmake_minimum_required (VERSION 2.6)
project (pmatch)

file (GLOB source src/*.cpp)
add_executable(pmatch ${source})
