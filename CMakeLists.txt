project(LBMgraphique)

cmake_minimum_required(VERSION 2.8.12)

set(CMAKE_PREFIX_PATH /opt/homebrew/opt/qt5/lib/cmake/)

find_package(Qt5Widgets REQUIRED)
find_package(Eigen3 REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
set(CMAKE_BUILD_TYPE Release)

set(SOURCES main.cpp grid.cpp)
    
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)


set(HEADERS grid.h lbm.hpp Debug.hpp LBMController.hpp settings.hpp)

add_executable(LBMgraphique
  ${SOURCES}
  ${HEADERS}
  )

target_include_directories(LBMgraphique
  PRIVATE
  ${QT_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  )

target_link_libraries( 
  LBMgraphique
  ${CMAKE_THREAD_LIBS_INIT}
  ${Qt5Widgets_LIBRARIES}
  )
