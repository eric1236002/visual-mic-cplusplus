cmake_minimum_required(VERSION 3.15)
project(VisualMicrophone VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Aggressive optimization flags for FFT performance
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native -ffast-math -funroll-loops")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Add local cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
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

# Link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(visual_microphone PUBLIC OpenMP::OpenMP_CXX)
endif()

# Installation
install(TARGETS visual_microphone DESTINATION bin)

