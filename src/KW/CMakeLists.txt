#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# Define this directory's FORTRAN flags:
if(WITH_ROOT)
  set(${PROJECT_NAME}_KW_FORTRAN_FLAGS "-g -mcmodel=large -fbacktrace -cpp -ffpe-trap=invalid,zero -Wall -fno-automatic -fcheck=bounds -O2 -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments -D__ROOT__" CACHE STRING "${PROJECT_NAME} KW FORTRAN FLAGS")
else()
  set(${PROJECT_NAME}_KW_FORTRAN_FLAGS "-g -mcmodel=large -fbacktrace -cpp -ffpe-trap=invalid,zero -Wall -fno-automatic -fcheck=bounds -O2 -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments -D__BS__" CACHE STRING "${PROJECT_NAME} KW FORTRAN FLAGS")
endif()

list(APPEND FORTRAN_FLAGS "${${PROJECT_NAME}_KW_FORTRAN_FLAGS}")
message(STATUS "fortran version ${CMAKE_CXX_COMPILER_VERSION}")
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0.0")
  list(APPEND FORTRAN_FLAGS "-fallow-argument-mismatch")
endif()
list(JOIN FORTRAN_FLAGS " " ${PROJECT_NAME}_KW_FORTRAN_FLAGS)

# Enable this directory's FORTRAN flags:
set(CMAKE_Fortran_FLAGS "${${PROJECT_NAME}_KW_FORTRAN_FLAGS}")
message(STATUS "fortran flag ${CMAKE_Fortran_FLAGS}")

# Define this directory's CXX flags:
set(${PROJECT_NAME}_KW_CXX_FLAGS "-O2 -g -Wall -fPIC -pthread -m64" CACHE STRING "${PROJECT_NAME} KW CXX FLAGS")
# Enable this directory's CXX flags:
set(CMAKE_CXX_FLAGS "${${PROJECT_NAME}_KW_CXX_FLAGS}")


# Do the local stuff:
set(KW_SRC
  bas.f
  bas.cpp
  bor.f
  conf.f
  d.f
  dky.f
  ems.f
  fra.f
  freeze.f
  graph2.f
  graph.f
  h.f
  hg.f
  hnb.f
  ho.h
  ico.f
  ids.f
  int.f
  jpo.f
  key.f
  lea.f  
  ../modif.cpp
  omg.f
  rsh.f
  sem.f
  tab.f
  tim.f
  u.f
  uc.f
  uti.f
  xan.f
  xht.f
)

if(WITH_ROOT)
  list(APPEND KW_SRC 
    eos.f
    hyt.f
    out.f)
endif()

add_library(KW OBJECT ${KW_SRC})
target_include_directories(
    KW 
    PRIVATE ../../src/TP/ 
    PRIVATE ../../src/KWb/ 
    PRIVATE ../../srcext/UR/)

if(WITH_ROOT)
  target_include_directories(
    KW
    PRIVATE ${ROOT_INCLUDE_DIR})
endif()

# Target XeposC
set(XEPOSC_SRC
  c.f
  d.f
)

add_library(XEPOSC OBJECT ${XEPOSC_SRC})
