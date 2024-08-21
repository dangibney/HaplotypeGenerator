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

# Utility rule file for benchmark.

# Include any custom commands dependencies for this target.
include tests/CMakeFiles/benchmark.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/benchmark.dir/progress.make

tests/CMakeFiles/benchmark:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/home/dan/Downloads/igraph-0.10.12/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Running benchmarks..."
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && true

benchmark: tests/CMakeFiles/benchmark
benchmark: tests/CMakeFiles/benchmark.dir/build.make
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: graphicality"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_graphicality
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_average_path_length_unweighted"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_average_path_length_unweighted
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_betweenness"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_betweenness
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_betweenness_weighted"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_betweenness_weighted
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_cliques"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_cliques
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_closeness_weighted"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_closeness_weighted
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_coloring"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_coloring
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_decompose"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_decompose
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_degree"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_degree
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_degree_sequence_game"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_degree_sequence_game
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_distances"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_distances
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_ecc"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_ecc
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_induced_subgraph_edges"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_induced_subgraph_edges
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_layout_umap"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_layout_umap
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_matrix_transpose"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_matrix_transpose
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_maximal_cliques"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_maximal_cliques
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_neighborhood"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_neighborhood
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_pagerank"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_pagerank
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_pagerank_weighted"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_pagerank_weighted
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_power_law_fit"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_power_law_fit
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_qsort"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_qsort
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_random_walk"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_random_walk
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_strength"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_strength
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_transitivity"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_transitivity
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_tree_game"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_tree_game
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_vertex_connectivity"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_vertex_connectivity
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: igraph_voronoi"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_igraph_voronoi
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold "Running benchmark: inc_vs_adj"
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && ./benchmark_inc_vs_adj
.PHONY : benchmark

# Rule to build all files generated by this target.
tests/CMakeFiles/benchmark.dir/build: benchmark
.PHONY : tests/CMakeFiles/benchmark.dir/build

tests/CMakeFiles/benchmark.dir/clean:
	cd /home/dan/Downloads/igraph-0.10.12/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/benchmark.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/benchmark.dir/clean

tests/CMakeFiles/benchmark.dir/depend:
	cd /home/dan/Downloads/igraph-0.10.12/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dan/Downloads/igraph-0.10.12 /home/dan/Downloads/igraph-0.10.12/tests /home/dan/Downloads/igraph-0.10.12/build /home/dan/Downloads/igraph-0.10.12/build/tests /home/dan/Downloads/igraph-0.10.12/build/tests/CMakeFiles/benchmark.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : tests/CMakeFiles/benchmark.dir/depend

