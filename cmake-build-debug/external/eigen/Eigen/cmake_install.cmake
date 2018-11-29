# Install script for directory: /home/bq2139/Documents/ogl/external/eigen/Eigen

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/Cholesky"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/CholmodSupport"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/Core"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/Dense"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/Eigen"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/Eigenvalues"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/Geometry"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/Householder"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/IterativeLinearSolvers"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/Jacobi"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/KLUSupport"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/LU"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/MetisSupport"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/OrderingMethods"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/PaStiXSupport"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/PardisoSupport"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/QR"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/QtAlignedMalloc"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/SPQRSupport"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/SVD"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/Sparse"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/SparseCholesky"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/SparseCore"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/SparseLU"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/SparseQR"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/StdDeque"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/StdList"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/StdVector"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/SuperLUSupport"
    "/home/bq2139/Documents/ogl/external/eigen/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/home/bq2139/Documents/ogl/external/eigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

