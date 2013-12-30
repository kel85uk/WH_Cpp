
include_directories( ../src )
add_executable( WH_test ../src/boundary_interior_tau_wh.cc ../src/double_compares.cpp meshfree_main.cc )
