cmake_minimum_required(VERSION 3.15)
project(VisualMicrophone VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Add local cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Find MPI package
find_package(MPI REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${MPI_CXX_INCLUDE_DIRS}
)

# Source files
set(SOURCES
    src/main.cpp
    src/steerable_pyramid.cpp
    src/signal_processing.cpp
    src/sound_from_video.cpp
    src/utils.cpp
    src/math_ops.cpp
    src/frame_loader.cpp
)

# Header files
set(HEADERS
    include/steerable_pyramid.h
    include/signal_processing.h
    include/sound_from_video.h
    include/utils.h
    include/math_ops.h
    include/frame_loader.h
)

# Create executable
add_executable(visual_microphone ${SOURCES} ${HEADERS})

# Link MPI libraries
target_link_libraries(visual_microphone 
    PUBLIC 
    MPI::MPI_CXX
)

# Add MPI compile flags
target_compile_options(visual_microphone 
    PRIVATE 
    ${MPI_CXX_COMPILE_FLAGS}
)

# Add MPI link flags
set_target_properties(visual_microphone 
    PROPERTIES 
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}"
)

# Installation
install(TARGETS visual_microphone DESTINATION bin)

