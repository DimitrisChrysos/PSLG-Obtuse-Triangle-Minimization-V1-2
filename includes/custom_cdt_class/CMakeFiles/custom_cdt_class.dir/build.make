# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/dimitris/Documents/project emiris/project2"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/dimitris/Documents/project emiris/project2"

# Include any dependencies generated for this target.
include includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/compiler_depend.make

# Include the progress variables for this target.
include includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/progress.make

# Include the compile flags for this target's objects.
include includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/flags.make

includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.o: includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/flags.make
includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.o: includes/custom_cdt_class/custom_cdt_class.cpp
includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.o: includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/dimitris/Documents/project emiris/project2/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.o"
	cd "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.o -MF CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.o.d -o CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.o -c "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class/custom_cdt_class.cpp"

includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.i"
	cd "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class/custom_cdt_class.cpp" > CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.i

includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.s"
	cd "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class/custom_cdt_class.cpp" -o CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.s

# Object files for target custom_cdt_class
custom_cdt_class_OBJECTS = \
"CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.o"

# External object files for target custom_cdt_class
custom_cdt_class_EXTERNAL_OBJECTS =

includes/custom_cdt_class/libcustom_cdt_class.a: includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/custom_cdt_class.cpp.o
includes/custom_cdt_class/libcustom_cdt_class.a: includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/build.make
includes/custom_cdt_class/libcustom_cdt_class.a: includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/dimitris/Documents/project emiris/project2/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libcustom_cdt_class.a"
	cd "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class" && $(CMAKE_COMMAND) -P CMakeFiles/custom_cdt_class.dir/cmake_clean_target.cmake
	cd "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/custom_cdt_class.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/build: includes/custom_cdt_class/libcustom_cdt_class.a
.PHONY : includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/build

includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/clean:
	cd "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class" && $(CMAKE_COMMAND) -P CMakeFiles/custom_cdt_class.dir/cmake_clean.cmake
.PHONY : includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/clean

includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/depend:
	cd "/home/dimitris/Documents/project emiris/project2" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/dimitris/Documents/project emiris/project2" "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class" "/home/dimitris/Documents/project emiris/project2" "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class" "/home/dimitris/Documents/project emiris/project2/includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : includes/custom_cdt_class/CMakeFiles/custom_cdt_class.dir/depend

