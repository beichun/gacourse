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
include external/eigen/blas/testing/CMakeFiles/cblat3.dir/depend.make

# Include the progress variables for this target.
include external/eigen/blas/testing/CMakeFiles/cblat3.dir/progress.make

# Include the compile flags for this target's objects.
include external/eigen/blas/testing/CMakeFiles/cblat3.dir/flags.make

external/eigen/blas/testing/CMakeFiles/cblat3.dir/cblat3.f.o: external/eigen/blas/testing/CMakeFiles/cblat3.dir/flags.make
external/eigen/blas/testing/CMakeFiles/cblat3.dir/cblat3.f.o: ../external/eigen/blas/testing/cblat3.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bq2139/Documents/ogl/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object external/eigen/blas/testing/CMakeFiles/cblat3.dir/cblat3.f.o"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/bq2139/Documents/ogl/external/eigen/blas/testing/cblat3.f -o CMakeFiles/cblat3.dir/cblat3.f.o

external/eigen/blas/testing/CMakeFiles/cblat3.dir/cblat3.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/cblat3.dir/cblat3.f.i"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/bq2139/Documents/ogl/external/eigen/blas/testing/cblat3.f > CMakeFiles/cblat3.dir/cblat3.f.i

external/eigen/blas/testing/CMakeFiles/cblat3.dir/cblat3.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/cblat3.dir/cblat3.f.s"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/bq2139/Documents/ogl/external/eigen/blas/testing/cblat3.f -o CMakeFiles/cblat3.dir/cblat3.f.s

# Object files for target cblat3
cblat3_OBJECTS = \
"CMakeFiles/cblat3.dir/cblat3.f.o"

# External object files for target cblat3
cblat3_EXTERNAL_OBJECTS =

external/eigen/blas/testing/cblat3: external/eigen/blas/testing/CMakeFiles/cblat3.dir/cblat3.f.o
external/eigen/blas/testing/cblat3: external/eigen/blas/testing/CMakeFiles/cblat3.dir/build.make
external/eigen/blas/testing/cblat3: external/eigen/blas/libeigen_blas.so
external/eigen/blas/testing/cblat3: external/eigen/blas/testing/CMakeFiles/cblat3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bq2139/Documents/ogl/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable cblat3"
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/blas/testing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cblat3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/eigen/blas/testing/CMakeFiles/cblat3.dir/build: external/eigen/blas/testing/cblat3

.PHONY : external/eigen/blas/testing/CMakeFiles/cblat3.dir/build

external/eigen/blas/testing/CMakeFiles/cblat3.dir/clean:
	cd /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/blas/testing && $(CMAKE_COMMAND) -P CMakeFiles/cblat3.dir/cmake_clean.cmake
.PHONY : external/eigen/blas/testing/CMakeFiles/cblat3.dir/clean

external/eigen/blas/testing/CMakeFiles/cblat3.dir/depend:
	cd /home/bq2139/Documents/ogl/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bq2139/Documents/ogl /home/bq2139/Documents/ogl/external/eigen/blas/testing /home/bq2139/Documents/ogl/cmake-build-debug /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/blas/testing /home/bq2139/Documents/ogl/cmake-build-debug/external/eigen/blas/testing/CMakeFiles/cblat3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/eigen/blas/testing/CMakeFiles/cblat3.dir/depend

