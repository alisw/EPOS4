#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

if(WITH_XA)
  set(XA_SRC 
    pack9.f
    z1.f
    z2.f 
    z3.cpp
    z4.c
    )
  # Define this directory's flags:
  set(${PROJECT_NAME}_XA_FORTRAN_FLAGS "-g -mcmodel=large -cpp -fbacktrace -ffpe-trap=invalid,zero  -Wall -fno-automatic -fcheck=bounds -O2 -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments" CACHE STRING "${PROJECT_NAME} XCTEQ FORTRAN FLAGS")
else()
  # VO is unused... link with a dummy interface
  set(XA_SRC
    dummyInterface.f
    dummyInterface.c
    )
  # Define this directory's flags:
  set(${PROJECT_NAME}_XA_FORTRAN_FLAGS "-g -mcmodel=large -cpp -fbacktrace -ffpe-trap=invalid,zero  -Wall -fno-automatic -fcheck=bounds -O2 -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments -Wno-unused-dummy-argument" CACHE STRING "${PROJECT_NAME} XCTEQ FORTRAN FLAGS")
endif()

# Enable this directory's flags:
set(CMAKE_Fortran_FLAGS "${${PROJECT_NAME}_XA_FORTRAN_FLAGS}")

# Define this directory's flags:
set(${PROJECT_NAME}_XA_CXX_FLAGS "-O2 -g -mcmodel=large -Wall -fPIC -pthread -m64" CACHE STRING "${PROJECT_NAME} XA CXX FLAGS")
# Enable this directory's flags:
set(CMAKE_CXX_FLAGS "${${PROJECT_NAME}_XA_CXX_FLAGS}")

# Define the project's flags:
set(${PROJECT_NAME}_XA_C_FLAGS "-g -mcmodel=large" CACHE STRING "${PROJECT_NAME} C FLAGS")
# Enable the project's flags:
set(CMAKE_C_FLAGS "${${PROJECT_NAME}_XA_C_FLAGS}")

#add_library(XA SHARED ${XA_SRC})
## set_target_properties(bg PROPERTIES LINKER_LANGUAGE FORTRAN)
#target_include_directories(XA PRIVATE ../KW/ ../VO/)
#install(TARGETS XA DESTINATION lib)

add_library(XA OBJECT ${XA_SRC})
target_include_directories(XA PRIVATE ../../src/KW/ ../../srcext/VO/)
