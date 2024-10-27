# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

# Define project name
project( polyg )


# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt5)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED COMPONENTS system filesystem)

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local package


# Creating entries for target: polyg
# ############################

add_executable( polyg  ex.cpp )

# Link the executable to CGAL and third-party libraries
target_link_libraries(polyg PRIVATE CGAL::CGAL ${EXTRA_LIBS})

# Link Qt5 libraries only if CGAL_Qt5 is found
if (CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER) # Enable viewer functionality
  target_link_libraries(polyg PRIVATE CGAL::CGAL_Qt5)
endif()