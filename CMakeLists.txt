cmake_minimum_required (VERSION 2.6)
project (Simulator)

list(APPEND CMAKE_CXX_FLAGS "-O3 -pthread")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Include library headers into the project
include_directories(lib/Matrix)

# Add the library source to the build
add_subdirectory(lib/Matrix)

# Include headers into the project
include_directories(include)

# Add the executable
add_executable(Simulator src/Driver.cpp src/Node.cpp include/Node.hpp src/Packet.cpp include/Packet.hpp src/Simulator.cpp include/Simulator.hpp src/NetGen.cpp include/Queue.hpp)

# Installation location
install(TARGETS Simulator DESTINATION bin)

# Link the library
target_link_libraries(Simulator Matrix)
