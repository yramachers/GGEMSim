# - Basic CMake setup
# Check version meets your requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.20)
project(THGEMSig)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# enable threading
list(APPEND CMAKE_CXX_FLAGS "-pthread -std=c++20 ${CMAKE_CXX_FLAGS}")

#Require ROOT, initially try finding previously installed root
find_package(ROOT REQUIRED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

# START of wirepulses part
include_directories(include ${ROOT_INCLUDES})

# Build a dynamic library from our sources
add_library(thgem SHARED 
  include/thread_pool.hpp
  include/fields.hh 
  include/geomodel.hh
  include/transport.hh 
  include/electrode.hh 
  src/thread_pool.cpp
  src/fields.cpp 
  src/geomodel.cpp 
  src/transport.cpp 
  src/electrode.cpp )
target_link_libraries(thgem ${ROOT_LIBRARIES})

#Executables
add_executable(signal_calculator.exe src/signal_calculator.cpp)
target_link_libraries(signal_calculator.exe ${ROOT_LIBRARIES} thgem ROOT::Geom)

#add_executable(signal_scan.exe src/signal_scan.cpp)
#target_link_libraries(signal_scan.exe ${ROOT_LIBRARIES} thgem ROOT::Geom)

#add_executable(multi_scan.exe src/multi_scan.cpp)
#target_link_libraries(multi_scan.exe ${ROOT_LIBRARIES} thgem ROOT::Geom)

#add_executable(testCS.exe utils/testCS.cpp)
#target_link_libraries(testCS.exe ${ROOT_LIBRARIES})
