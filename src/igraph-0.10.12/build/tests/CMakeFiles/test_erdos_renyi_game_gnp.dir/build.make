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
include tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/flags.make

tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.o: tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/flags.make
tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.o: /home/dan/Downloads/igraph-0.10.12/tests/unit/erdos_renyi_game_gnp.c
tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.o: tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dan/Downloads/igraph-0.10.12/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.o"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.o -MF CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.o.d -o CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.o -c /home/dan/Downloads/igraph-0.10.12/tests/unit/erdos_renyi_game_gnp.c

tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.i"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/dan/Downloads/igraph-0.10.12/tests/unit/erdos_renyi_game_gnp.c > CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.i

tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.s"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/dan/Downloads/igraph-0.10.12/tests/unit/erdos_renyi_game_gnp.c -o CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.s

# Object files for target test_erdos_renyi_game_gnp
test_erdos_renyi_game_gnp_OBJECTS = \
"CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.o"

# External object files for target test_erdos_renyi_game_gnp
test_erdos_renyi_game_gnp_EXTERNAL_OBJECTS = \
"/home/dan/Downloads/igraph-0.10.12/build/tests/CMakeFiles/test_utilities.dir/unit/test_utilities.c.o"

tests/test_erdos_renyi_game_gnp: tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/unit/erdos_renyi_game_gnp.c.o
tests/test_erdos_renyi_game_gnp: tests/CMakeFiles/test_utilities.dir/unit/test_utilities.c.o
tests/test_erdos_renyi_game_gnp: tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/build.make
tests/test_erdos_renyi_game_gnp: src/libigraph.a
tests/test_erdos_renyi_game_gnp: /usr/lib/x86_64-linux-gnu/libm.so
tests/test_erdos_renyi_game_gnp: /usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so
tests/test_erdos_renyi_game_gnp: /usr/lib/x86_64-linux-gnu/libpthread.a
tests/test_erdos_renyi_game_gnp: tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/dan/Downloads/igraph-0.10.12/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_erdos_renyi_game_gnp"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_erdos_renyi_game_gnp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/build: tests/test_erdos_renyi_game_gnp
.PHONY : tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/build

tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/clean:
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/test_erdos_renyi_game_gnp.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/clean

tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/depend:
	cd /home/dan/Downloads/igraph-0.10.12/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dan/Downloads/igraph-0.10.12 /home/dan/Downloads/igraph-0.10.12/tests /home/dan/Downloads/igraph-0.10.12/build /home/dan/Downloads/igraph-0.10.12/build/tests /home/dan/Downloads/igraph-0.10.12/build/tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : tests/CMakeFiles/test_erdos_renyi_game_gnp.dir/depend

