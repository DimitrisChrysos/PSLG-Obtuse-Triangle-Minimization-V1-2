# Create the utils library
add_library(custom_cdt_class custom_cdt_class.cpp)

# Include the current directory for headers
target_include_directories(custom_cdt_class PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})