# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name and supported languages
project(DNA_Hidden_Repeat_Detection LANGUAGES CXX)

# ✅ Set C++ standard to C++17 (or higher if needed)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ✅ Set output directories for binary files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ✅ Include header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ✅ Find all source files
file(GLOB_RECURSE SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/*.h
)

# ✅ Exclude CMake-generated files from the build
list(FILTER SOURCES EXCLUDE REGEX ".*CMakeFiles.*")
list(FILTER SOURCES EXCLUDE REGEX ".*CompilerIdCXX.*")

# ✅ Ensure that sources are found, otherwise throw an error
if (NOT SOURCES)
    message(FATAL_ERROR "No source files found! Check paths or file names.")
endif()

# ✅ Define the executable name
add_executable(dna-hidden-repeat-detector ${SOURCES})

# ✅ Add debug/release flags
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Building in Debug Mode")
    target_compile_definitions(dna-hidden-repeat-detector PRIVATE DEBUG)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    message("Building in Release Mode")
    target_compile_definitions(dna-hidden-repeat-detector PRIVATE NDEBUG)
endif()

# ✅ Enable strict warnings
if (MSVC)
    target_compile_options(dna-hidden-repeat-detector PRIVATE /W4)
else()
    target_compile_options(dna-hidden-repeat-detector PRIVATE -Wall -Wextra -pedantic)
endif()

# ✅ Set optimization flags for Release builds
if(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

# ✅ Windows-specific configuration
if(WIN32)
    message(STATUS "Building for Windows")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    target_compile_definitions(dna-hidden-repeat-detector PRIVATE _WIN32)
endif()

# ✅ Install target for 'make install'
install(TARGETS dna-hidden-repeat-detector
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# ✅ Clean command (optional)
add_custom_target(clean-all
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
