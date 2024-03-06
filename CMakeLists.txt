cmake_minimum_required(VERSION 3.25)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(app)

set(CMAKE_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX14.sdk")

set(Matlab_ROOT_DIR "/Applications/MATLAB_R2023b.app")

# Find MATLAB
find_package(Matlab REQUIRED MX_LIBRARY MAT_LIBRARY ENG_LIBRARY)

find_package(Eigen3 REQUIRED)
# Include directories
target_include_directories(app PRIVATE ${Matlab_INCLUDE_DIRS})
target_include_directories(app PRIVATE ${EIGEN3_INCLUDE_DIR})


file(GLOB OBJECT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/*.h)

set(SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/*.cpp
        ${OBJECT_HEADERS}
)

add_executable(app ${SOURCE_FILES})
# Define data directories
set(DATA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Data")
set(DATA_DEST_DIR "${CMAKE_BINARY_DIR}/Data")

# Custom target for copying data
add_custom_target(data
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${DATA_SRC_DIR} ${DATA_DEST_DIR}
    )
# Create executable

add_dependencies(app data)

# Link MATLAB libraries
target_link_libraries(app ${Matlab_LIBRARIES})














