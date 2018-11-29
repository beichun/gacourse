# Install script for directory: /home/bq2139/Documents/ogl/external/assimp-3.0.1270/code

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibassimp3.0-r1270-OGLtutsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/bq2139/Documents/ogl/cmake-build-debug/external/assimp-3.0.1270/code/libassimp.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/assimp" TYPE FILE FILES
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/anim.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/ai_assert.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/camera.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/color4.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/color4.inl"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/config.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/defs.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/cfileio.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/light.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/material.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/material.inl"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/matrix3x3.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/matrix3x3.inl"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/matrix4x4.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/matrix4x4.inl"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/mesh.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/postprocess.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/quaternion.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/quaternion.inl"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/scene.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/texture.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/types.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/vector2.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/vector2.inl"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/vector3.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/vector3.inl"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/version.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/cimport.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/importerdesc.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/Importer.hpp"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/DefaultLogger.hpp"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/ProgressHandler.hpp"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/IOStream.hpp"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/IOSystem.hpp"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/Logger.hpp"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/LogStream.hpp"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/NullLogger.hpp"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/cexport.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/Exporter.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/assimp/Compiler" TYPE FILE FILES
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/Compiler/pushpack1.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/../include/assimp/Compiler/poppack1.h"
    "/home/bq2139/Documents/ogl/external/assimp-3.0.1270/code/pstdint.h"
    )
endif()

