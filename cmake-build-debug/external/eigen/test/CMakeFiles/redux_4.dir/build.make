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
include external/eigen/test/CMakeFiles/redux_4.dir/depend.make

# Include the progress variables for this target.
include external/eigen/test/CMakeFiles/redux_4.dir/progress.make

# Include the compile flags for this target's objects.
include external/eigen/test/CMakeFiles/redux_4.dir/flags.make

external/eigen/test/CMakeFiles/redux_4.dir/redux.cpp.o: external/eigen/test/CMakeFiles/redux_4.dir/flags.make
external/eigen/test/CMakeFiles/redux_4.dir/redux.cpp.o: ../external/eigen/test/redux.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bq2139/Documents/ogl/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/eigen/test/CMakeFiles/redux_4.dir/redux.cpp.o"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/redux_4.dir/redux.cpp.o -c /home/bq2139/Documents/ogl/external/eigen/test/redux.cpp

external/eigen/test/CMakeFiles/redux_4.dir/redux.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/redux_4.dir/redux.cpp.i"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/bq2139/Documents/ogl/external/eigen/test/redux.cpp > CMakeFiles/redux_4.dir/redux.cpp.i

external/eigen/test/CMakeFiles/redux_4.dir/redux.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/redux_4.dir/redux.cpp.s"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/bq2139/Documents/ogl/external/eigen/test/redux.cpp -o CMakeFiles/redux_4.dir/redux.cpp.s

# Object files for target redux_4
redux_4_OBJECTS = \
"CMakeFiles/redux_4.dir/redux.cpp.o"

# External object files for target redux_4
redux_4_EXTERNAL_OBJECTS =

external/eigen/test/redux_4: external/eigen/test/CMakeFiles/redux_4.dir/redux.cpp.o
external/eigen/test/redux_4: external/eigen/test/CMakeFiles/redux_4.dir/build.make
external/eigen/test/redux_4: external/eigen/test/CMakeFiles/redux_4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bq2139/Documents/ogl/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable redux_4"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/redux_4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/eigen/test/CMakeFiles/redux_4.dir/build: external/eigen/test/redux_4

.PHONY : external/eigen/test/CMakeFiles/redux_4.dir/build

external/eigen/test/CMakeFiles/redux_4.dir/clean:
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/test && $(CMAKE_COMMAND) -P CMakeFiles/redux_4.dir/cmake_clean.cmake
.PHONY : external/eigen/test/CMakeFiles/redux_4.dir/clean

external/eigen/test/CMakeFiles/redux_4.dir/depend:
	cd /home/bq2139/Documents/ogl/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bq2139/Documents/ogl /home/bq2139/Documents/ogl/external/eigen/test /home/bq2139/Documents/ogl/cmake-build-debug /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/test /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/test/CMakeFiles/redux_4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/eigen/test/CMakeFiles/redux_4.dir/depend

