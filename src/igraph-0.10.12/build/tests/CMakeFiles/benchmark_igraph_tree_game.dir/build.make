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

# Include any dependencies generated for this target.
include tests/CMakeFiles/benchmark_igraph_tree_game.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/CMakeFiles/benchmark_igraph_tree_game.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/benchmark_igraph_tree_game.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/benchmark_igraph_tree_game.dir/flags.make

tests/CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.o: tests/CMakeFiles/benchmark_igraph_tree_game.dir/flags.make
tests/CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.o: /home/dan/Downloads/igraph-0.10.12/tests/benchmarks/igraph_tree_game.c
tests/CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.o: tests/CMakeFiles/benchmark_igraph_tree_game.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dan/Downloads/igraph-0.10.12/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object tests/CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.o"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT tests/CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.o -MF CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.o.d -o CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.o -c /home/dan/Downloads/igraph-0.10.12/tests/benchmarks/igraph_tree_game.c

tests/CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.i"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/dan/Downloads/igraph-0.10.12/tests/benchmarks/igraph_tree_game.c > CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.i

tests/CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.s"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/dan/Downloads/igraph-0.10.12/tests/benchmarks/igraph_tree_game.c -o CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.s

# Object files for target benchmark_igraph_tree_game
benchmark_igraph_tree_game_OBJECTS = \
"CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.o"

# External object files for target benchmark_igraph_tree_game
benchmark_igraph_tree_game_EXTERNAL_OBJECTS =

tests/benchmark_igraph_tree_game: tests/CMakeFiles/benchmark_igraph_tree_game.dir/benchmarks/igraph_tree_game.c.o
tests/benchmark_igraph_tree_game: tests/CMakeFiles/benchmark_igraph_tree_game.dir/build.make
tests/benchmark_igraph_tree_game: src/libigraph.a
tests/benchmark_igraph_tree_game: /usr/lib/x86_64-linux-gnu/libm.so
tests/benchmark_igraph_tree_game: /usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so
tests/benchmark_igraph_tree_game: /usr/lib/x86_64-linux-gnu/libpthread.a
tests/benchmark_igraph_tree_game: tests/CMakeFiles/benchmark_igraph_tree_game.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/dan/Downloads/igraph-0.10.12/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable benchmark_igraph_tree_game"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/benchmark_igraph_tree_game.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/benchmark_igraph_tree_game.dir/build: tests/benchmark_igraph_tree_game
.PHONY : tests/CMakeFiles/benchmark_igraph_tree_game.dir/build

tests/CMakeFiles/benchmark_igraph_tree_game.dir/clean:
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/benchmark_igraph_tree_game.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/benchmark_igraph_tree_game.dir/clean

tests/CMakeFiles/benchmark_igraph_tree_game.dir/depend:
	cd /home/dan/Downloads/igraph-0.10.12/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dan/Downloads/igraph-0.10.12 /home/dan/Downloads/igraph-0.10.12/tests /home/dan/Downloads/igraph-0.10.12/build /home/dan/Downloads/igraph-0.10.12/build/tests /home/dan/Downloads/igraph-0.10.12/build/tests/CMakeFiles/benchmark_igraph_tree_game.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : tests/CMakeFiles/benchmark_igraph_tree_game.dir/depend

