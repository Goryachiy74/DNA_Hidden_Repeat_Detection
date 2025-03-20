# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name and supported languages
project(DNA_Hidden_Repeat_Detection LANGUAGES CXX)

# Set C++ standard to C++17 (or C++20 if preferred)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set output directories for binary files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include header files
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find all source files automatically
file(GLOB_RECURSE SOURCES
    src/*.cpp
    include/*.h
)

# Define the executable name
add_executable(dna-hidden-repeat-detector ${SOURCES})

# Add debug/release flags
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Building in Debug Mode")
    target_compile_definitions(dna-hidden-repeat-detector PRIVATE DEBUG)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    message("Building in Release Mode")
    target_compile_definitions(dna-hidden-repeat-detector PRIVATE NDEBUG)
endif()

# Enable compiler warnings for stricter code checking
if (MSVC)
    target_compile_options(dna-hidden-repeat-detector PRIVATE /W4)
else()
    target_compile_options(dna-hidden-repeat-detector PRIVATE -Wall -Wextra -pedantic)
endif()

# Link required libraries (if any)
# target_link_libraries(dna-hidden-repeat-detector PRIVATE <library_name>)

# Optional: Set optimization flags for Release builds
if(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

# Install target for 'make install' (optional)
install(TARGETS dna-hidden-repeat-detector
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Add 'make clean' command (optional)
add_custom_target(clean-all
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
