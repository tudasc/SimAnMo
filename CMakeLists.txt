cmake_minimum_required(VERSION 3.13)

# set the project name and version
project(SimAnMo VERSION 1.5)

# Add new build types
message("* Adding build types...")
SET(CMAKE_CXX_FLAGS_ULTRAFAST
    "-ofast -DNDEBUG"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE 
)

#set(CMAKE_BUILD_TYPE "ULTRAFAST1")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(SimulatedAnnealingExtraP/LibraryInterface.h.in SimulatedAnnealingExtraP/VersionInfo.h)

#list(APPEND EXTRA_LIBS AlgLib)
list(APPEND EXTRA_LIBS SimAnMoLib)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/SimulatedAnnealingExtraP ${PROJECT_SOURCE_DIR}/SimulatedAnnealingExtraP/alglib")

# add the libraries
add_subdirectory(SimulatedAnnealingExtraP)
#add_subdirectory(SimulatedAnnealingExtraP/alglib)

# add the executable
add_executable(SimAnMo "")
add_executable(SimAnMoLibraryTester "")

add_subdirectory(SimAnMoDriver)
add_subdirectory(SimAnMoTest)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(SimAnMo PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(SimAnMo PUBLIC 
	SimAnMoLib
#	AlgLib	
)

target_link_libraries(SimAnMoLibraryTester PUBLIC 
	SimAnMoLib
#	AlgLib	
)

message(STATUS "foo include dir: ${CMAKE_BUILD_TYPE}")

# add the binary tree to the search path for include files
target_include_directories(SimAnMo PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           ${PROJECT_SOURCE_DIR}/SimulatedAnnealingExtraP
 )
