# Create the steiner_methods library
add_library(steiner_methods steiner_methods.cpp)

# Include the current directory for headers
target_include_directories(steiner_methods PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})