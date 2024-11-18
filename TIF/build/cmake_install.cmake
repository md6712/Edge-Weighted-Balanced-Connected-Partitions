# Install script for directory: C:/Users/Morteza/Downloads/lemon-1.3.1/lemon-1.3.1

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files/LEMON")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "E:/OneDriveSkema/OneDrive - SKEMA Business School/Data/Research/1- OnGoing/Graph - Balanced Tree/Heavest Balanced Tree/TIF/build/cmake/LEMONConfig.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("E:/OneDriveSkema/OneDrive - SKEMA Business School/Data/Research/1- OnGoing/Graph - Balanced Tree/Heavest Balanced Tree/TIF/build/lemon/cmake_install.cmake")
  include("E:/OneDriveSkema/OneDrive - SKEMA Business School/Data/Research/1- OnGoing/Graph - Balanced Tree/Heavest Balanced Tree/TIF/build/contrib/cmake_install.cmake")
  include("E:/OneDriveSkema/OneDrive - SKEMA Business School/Data/Research/1- OnGoing/Graph - Balanced Tree/Heavest Balanced Tree/TIF/build/demo/cmake_install.cmake")
  include("E:/OneDriveSkema/OneDrive - SKEMA Business School/Data/Research/1- OnGoing/Graph - Balanced Tree/Heavest Balanced Tree/TIF/build/tools/cmake_install.cmake")
  include("E:/OneDriveSkema/OneDrive - SKEMA Business School/Data/Research/1- OnGoing/Graph - Balanced Tree/Heavest Balanced Tree/TIF/build/doc/cmake_install.cmake")
  include("E:/OneDriveSkema/OneDrive - SKEMA Business School/Data/Research/1- OnGoing/Graph - Balanced Tree/Heavest Balanced Tree/TIF/build/test/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
  file(WRITE "E:/OneDriveSkema/OneDrive - SKEMA Business School/Data/Research/1- OnGoing/Graph - Balanced Tree/Heavest Balanced Tree/TIF/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
