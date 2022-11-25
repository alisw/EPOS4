#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# find_package(ROOT)
# Define this directory's flags:
set(${PROJECT_NAME}_KWb_FORTRAN_FLAGS "-g -mcmodel=large -fPIC -cpp -fbacktrace -ffpe-trap=invalid,zero -Wall -fno-automatic -fcheck=bounds -O2 -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments" CACHE STRING "${PROJECT_NAME} KWb FORTRAN FLAGS")
# Enable this directory's flags:
set(CMAKE_Fortran_FLAGS "${${PROJECT_NAME}_KWb_FORTRAN_FLAGS}")

# Define this directory's flags:
if(WITH_ROOT)
# set(${PROJECT_NAME}_KWb_CXX_FLAGS "-std=c++11 -O2 -g -Wall -fPIC -pthread -m64 -D__ROOT__" CACHE STRING "${PROJECT_NAME} KWb CXX FLAGS")

list(APPEND CXX_FLAGS ${ROOT_CPPFLAGS})
list(APPEND CXX_FLAGS "-O2 -g -Wall -fPIC -D__ROOT__")
list(JOIN CXX_FLAGS " " ${PROJECT_NAME}_KWb_CXX_FLAGS)

else()
set(${PROJECT_NAME}_KWb_CXX_FLAGS "-std=c++11 -O2 -g -Wall -fPIC -pthread -m64" CACHE STRING "${PROJECT_NAME} KWb CXX FLAGS")
endif()
# Enable this directory's flags:
set(CMAKE_CXX_FLAGS "${${PROJECT_NAME}_KWb_CXX_FLAGS}")
message(STATUS "KWb CXX FLAGS: " ${CMAKE_CXX_FLAGS})


# Do the local stuff:
set(KWb_SRC 
  af.f 
  al.cpp
  time.c
  x_s.f
)

if(WITH_ROOT)
  list(APPEND KWb_SRC 
    etreef.cpp
    tre.cpp
    etree.cpp
    fill_hep.cpp)
endif()

add_library(KWb OBJECT ${KWb_SRC})
target_include_directories(KWb PRIVATE ../KW/ PRIVATE)

if(WITH_ROOT)
  target_include_directories(KWb PRIVATE ${ROOT_INCLUDE_DIR})
endif()

