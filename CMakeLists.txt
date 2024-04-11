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

if(APPLE)
    # set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    # set(OpenMP_CXX_FLAGS "-Xpreprocessor -O3 -march=native ")
    # set(OpenMP_CXX_FLAGS "-Xpreprocessor -O3 -march=native -flto -ffast-math -fno-rtti")
    # set(OpenMP_CXX_FLAGS "-Xpreprocessor -O3 -march=native -flto -fno-rtti -fopenmp -I/opt/homebrew/opt/libomp/include")

    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    # set(OpenMP_CXX_LIB_NAMES "omp")
    # set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/estimatedRobustPoints.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/estimatedViews.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/Helper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/estimatedPoints.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/Refinement.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation/diagnosis.cpp

    # List any other .cpp file here
)

add_executable(app ${SOURCE_FILES})


target_compile_options(app PRIVATE -O3 ${OpenMP_CXX_FLAGS})
target_include_directories(app PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/implementation
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${Matlab_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    /opt/homebrew/opt/libomp/include
) # Also add any subfolders of implementaton and include here


add_dependencies(app data)

target_link_libraries(app ${Matlab_LIBRARIES} )
if(OpenMP_CXX_FLAGS)
    target_link_libraries(app ${OpenMP_omp_LIBRARY})
endif()
# I feel like you should also do "target_link_libraries(app ${EIGEN3_LIBRARIES})" unless its header only?


# Create target data
set(DATA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Data")
set(DATA_DEST_DIR "${CMAKE_BINARY_DIR}/Data")

add_custom_target(data
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${DATA_SRC_DIR} ${DATA_DEST_DIR}
)


