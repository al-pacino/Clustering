cmake_minimum_required(VERSION 2.5)

project(PamMPI)

find_package(MPI REQUIRED)

include_directories( src_nothreads )
include_directories( ${MPI_INCLUDE_PATH} )

set(SOURCE
	src_nothreads/MpiSupport.cpp
	src_nothreads/main.cpp)

add_executable(pam ${SOURCE})

target_link_libraries(pam ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(pam PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(pam PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
