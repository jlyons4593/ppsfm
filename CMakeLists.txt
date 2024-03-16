# Set minimum versions for CMake and C++
cmake_minimum_required(VERSION 3.20...3.25)
set(CMAKE_CXX_STANDARD 17)

# Export compile commands for LSPs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set application root directories
set(CMAKE_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX14.sdk")
set(Matlab_ROOT_DIR "/Applications/MATLAB_R2023b.app")

project(app)

# Find packages
find_package(Matlab REQUIRED MX_LIBRARY MAT_LIBRARY ENG_LIBRARY)
find_package(Eigen3 REQUIRED)


# Create target app
set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/ppsfm_initialisation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/pipelineManager.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/EstimatedFundamentalMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/pyramidalVisibilityScore.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/pairAffinityCalculator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/factorCompletion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/estimatedRobustViews.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/estimatedViews.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/Helper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/Refinement.cpp

    # List any other .cpp file here
)

add_executable(app ${SOURCE_FILES})

target_include_directories(app PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${Matlab_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
) # Also add any subfolders of implementaton and include here

add_dependencies(app data)
target_link_libraries(app ${Matlab_LIBRARIES})
# I feel like you should also do "target_link_libraries(app ${EIGEN3_LIBRARIES})" unless its header only?


# Create target data
set(DATA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Data")
set(DATA_DEST_DIR "${CMAKE_BINARY_DIR}/Data")

add_custom_target(data
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${DATA_SRC_DIR} ${DATA_DEST_DIR}
)












