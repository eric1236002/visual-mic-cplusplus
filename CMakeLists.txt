cmake_minimum_required(VERSION 3.15)
project(VisualMicrophone VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add local cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

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


# Installation
install(TARGETS visual_microphone DESTINATION bin)

