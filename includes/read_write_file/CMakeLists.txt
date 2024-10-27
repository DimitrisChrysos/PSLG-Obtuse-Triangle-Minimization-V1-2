# Create the utils library
add_library(read_write_file read_write_file.cpp)

# Include the current directory for headers
target_include_directories(read_write_file PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})