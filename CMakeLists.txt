cmake_minimum_required( VERSION 2.6 ) 

project ( WH_asme )

set( CMAKE_CXX_COMPILER "mpicxx")
set( CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wconversion -Wshadow" )
set( CMAKE_BUILD_TYPE "Debug" )

add_subdirectory( src   )
add_subdirectory( include   )
add_subdirectory( exec )
