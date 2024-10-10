# Create the utils library
add_library(graphics graphics.cpp)

# Include the current directory for headers
target_include_directories(graphics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
