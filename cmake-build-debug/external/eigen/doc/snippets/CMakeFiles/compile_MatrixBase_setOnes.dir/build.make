# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/bq2139/Documents/clion-2018.2.6/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/bq2139/Documents/clion-2018.2.6/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/bq2139/Documents/ogl

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/bq2139/Documents/ogl/cmake-build-debug

# Include any dependencies generated for this target.
include external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/depend.make

# Include the progress variables for this target.
include external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/progress.make

# Include the compile flags for this target's objects.
include external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/flags.make

external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.o: external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/flags.make
external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.o: external/eigen/doc/snippets/compile_MatrixBase_setOnes.cpp
external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.o: ../external/eigen/doc/snippets/MatrixBase_setOnes.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bq2139/Documents/ogl/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.o"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.o -c /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets/compile_MatrixBase_setOnes.cpp

external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.i"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets/compile_MatrixBase_setOnes.cpp > CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.i

external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.s"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets/compile_MatrixBase_setOnes.cpp -o CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.s

# Object files for target compile_MatrixBase_setOnes
compile_MatrixBase_setOnes_OBJECTS = \
"CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.o"

# External object files for target compile_MatrixBase_setOnes
compile_MatrixBase_setOnes_EXTERNAL_OBJECTS =

external/eigen/doc/snippets/compile_MatrixBase_setOnes: external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/compile_MatrixBase_setOnes.cpp.o
external/eigen/doc/snippets/compile_MatrixBase_setOnes: external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/build.make
external/eigen/doc/snippets/compile_MatrixBase_setOnes: external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bq2139/Documents/ogl/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable compile_MatrixBase_setOnes"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_MatrixBase_setOnes.dir/link.txt --verbose=$(VERBOSE)
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets && ./compile_MatrixBase_setOnes >/home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets/MatrixBase_setOnes.out

# Rule to build all files generated by this target.
external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/build: external/eigen/doc/snippets/compile_MatrixBase_setOnes

.PHONY : external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/build

external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/clean:
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets && $(CMAKE_COMMAND) -P CMakeFiles/compile_MatrixBase_setOnes.dir/cmake_clean.cmake
.PHONY : external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/clean

external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/depend:
	cd /home/bq2139/Documents/ogl/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bq2139/Documents/ogl /home/bq2139/Documents/ogl/external/eigen/doc/snippets /home/bq2139/Documents/ogl/cmake-build-debug /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/eigen/doc/snippets/CMakeFiles/compile_MatrixBase_setOnes.dir/depend

