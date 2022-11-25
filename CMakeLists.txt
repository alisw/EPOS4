#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 2.8.12.2...3.17.5)

# This is your project statement. You should always list languages.
project(EPOS Fortran C CXX)


# Compilation option FULL or not FULL
set(COMPILE_OPTION $ENV{COP})

if(COMPILE_OPTION STREQUAL "FULL")
  set(WITH_ROOT True)
  set(WITH_FASTJET True)
  set(WITH_HYDRODYNAMIC True)
  set(WITH_EB False)
  set(WITH_VO True)
  set(WITH_CTEQ True)
  set(WITH_XA True)
  set(WITH_JT True)
elseif(COMPILE_OPTION STREQUAL "DEMO") 
  set(WITH_ROOT False)
  set(WITH_FASTJET False)
  set(WITH_HYDRODYNAMIC False)
  set(WITH_EB False)
  set(WITH_VO False)
  set(WITH_CTEQ False)
  set(WITH_XA False)
  set(WITH_JT False)
elseif(COMPILE_OPTION STREQUAL "BASIC") 
  set(WITH_ROOT True) 
  set(WITH_FASTJET True)
  set(WITH_HYDRODYNAMIC True)
  set(WITH_EB False)
  set(WITH_VO False)
  set(WITH_CTEQ False)
  set(WITH_XA False)
  set(WITH_JT False)
elseif(COMPILE_OPTION STREQUAL "PHSD") 
  set(WITH_ROOT True) 
  set(WITH_FASTJET True)
  set(WITH_HYDRODYNAMIC False)
  set(WITH_EB True)
  set(WITH_VO False)
  set(WITH_CTEQ False)
  set(WITH_XA False)
  set(WITH_JT False)
endif()

message(STATUS "Compilation with ${COMPILE_OPTION} option")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules" CACHE PATH "Module Path" FORCE)

if(WITH_ROOT)
  message(STATUS "Looking for ROOT")
  find_package(ROOT)
  message(STATUS "ROOT version ${ROOT_VERSION}")
  message(STATUS "ROOT libraries ${ROOT_LIBRARIES}")
endif()

if(WITH_FASTJET)
  message(STATUS "Looking for FastJet")
  set(FASTJET_ROOT_DIR $ENV{FASTSYS})
  find_package(FastJet)
  message(STATUS "FASTJET version ${FASTJET_VERSION}")
  message(STATUS "FASTJET libraries ${FASTJET_LIBRARIES}")
endif()

# Fortran module files are created in $LIBDIR/Mod                                          
set(LIBRARY_DIRECTORY $ENV{LIBDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${LIBRARY_DIRECTORY}Mod)

# Define the project's flags:
set(${PROJECT_NAME}_FORTRAN_FLAGS "-g -mcmodel=large -cpp -fbacktrace -ffpe-trap=invalid,zero -Wall -fno-automatic -fcheck=bounds -O2 -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments" CACHE STRING "${PROJECT_NAME} FORTRAN FLAGS")
# Enable this directory's flags:
set(CMAKE_Fortran_FLAGS "${${PROJECT_NAME}_FORTRAN_FLAGS}")

## Enter the subdirectories from src:
add_subdirectory(src/KW)
add_subdirectory(src/KWb)
add_subdirectory(src/KWb/HepMC)
add_subdirectory(src/TP)

## Enter the subdirectories from srcext:
add_subdirectory(srcext/CTEQ) 
add_subdirectory(srcext/UR)
add_subdirectory(srcext/VO)
add_subdirectory(srcext/XA)

if(WITH_JT)
  add_subdirectory(srcext/JT)
endif()

if(WITH_HYDRODYNAMIC)
  add_subdirectory(src/MS)
  add_subdirectory(srcext/YK)
endif()


####################################################################################
####################################################################################
# Xepos target
####################################################################################
####################################################################################
set(OBJECT_LIST
  $<TARGET_OBJECTS:KW> 
  $<TARGET_OBJECTS:KWb> 
  $<TARGET_OBJECTS:HEPMC> 
  $<TARGET_OBJECTS:TP> 
  $<TARGET_OBJECTS:UR> 
  $<TARGET_OBJECTS:CTEQ>	
  $<TARGET_OBJECTS:VO> 	
  $<TARGET_OBJECTS:XA> 
  )

if(WITH_JT)
  list(APPEND OBJECT_LIST
    $<TARGET_OBJECTS:JT> 
    )
endif()

if(WITH_HYDRODYNAMIC)
  list(APPEND OBJECT_LIST
    $<TARGET_OBJECTS:MS> 
    $<TARGET_OBJECTS:YK> 
    )
endif()

add_dependencies(XA VO)

add_executable(Xepos src/main.cpp ${OBJECT_LIST})
set_target_properties(Xepos PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(Xepos PROPERTIES LINK_FLAGS "-pthread -rdynamic -lgfortran")

set(LINK_LIBRARY_LIST "")
if(WITH_ROOT)
  list(APPEND LINK_LIBRARY_LIST
    ${ROOT_LIBRARIES}
    )
endif()
if(WITH_FASTJET)
  list(APPEND LINK_LIBRARY_LIST
    ${FASTJET_LIBRARIES}
    )
endif()

target_link_libraries(Xepos ${LINK_LIBRARY_LIST})

####################################################################################
####################################################################################
# XeposC target
####################################################################################
####################################################################################
add_executable(XeposC $<TARGET_OBJECTS:XEPOSC>)
set_target_properties(XeposC PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(XeposC PROPERTIES LINK_FLAGS "-pthread -rdynamic")
