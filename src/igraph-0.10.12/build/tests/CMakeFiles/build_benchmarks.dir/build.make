# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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
CMAKE_SOURCE_DIR = /home/dan/Downloads/igraph-0.10.12

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dan/Downloads/igraph-0.10.12/build

# Utility rule file for build_benchmarks.

# Include any custom commands dependencies for this target.
include tests/CMakeFiles/build_benchmarks.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/build_benchmarks.dir/progress.make

build_benchmarks: tests/CMakeFiles/build_benchmarks.dir/build.make
.PHONY : build_benchmarks

# Rule to build all files generated by this target.
tests/CMakeFiles/build_benchmarks.dir/build: build_benchmarks
.PHONY : tests/CMakeFiles/build_benchmarks.dir/build

tests/CMakeFiles/build_benchmarks.dir/clean:
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/build_benchmarks.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/build_benchmarks.dir/clean

tests/CMakeFiles/build_benchmarks.dir/depend:
	cd /home/dan/Downloads/igraph-0.10.12/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dan/Downloads/igraph-0.10.12 /home/dan/Downloads/igraph-0.10.12/tests /home/dan/Downloads/igraph-0.10.12/build /home/dan/Downloads/igraph-0.10.12/build/tests /home/dan/Downloads/igraph-0.10.12/build/tests/CMakeFiles/build_benchmarks.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : tests/CMakeFiles/build_benchmarks.dir/depend

