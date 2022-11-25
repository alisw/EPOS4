#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# find_package(FastJet)
# Define this directory's flags:
set(${PROJECT_NAME}_TP_FORTRAN_FLAGS "-g -mcmodel=large -cpp -fbacktrace -ffpe-trap=invalid,zero  -Wall -fno-automatic -fcheck=bounds -O2 -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments" CACHE STRING "${PROJECT_NAME} TP FORTRAN FLAGS")
# Enable this directory's flags:
set(CMAKE_Fortran_FLAGS "${${PROJECT_NAME}_TP_FORTRAN_FLAGS}")

# Define this directory's flags:
set(${PROJECT_NAME}_TP_CXX_FLAGS "-g -O2 -Wall -fPIC" CACHE STRING "${PROJECT_NAME} JT CXX FLAGS")
# Enable this directory's flags:
set(CMAKE_CXX_FLAGS "${${PROJECT_NAME}_TP_CXX_FLAGS}")

# Do the local stuff:
set(TP_SRC 
  gen.f
  models.f  
  par.f 
  pdf.f
  ptb.f
  qsh.f
  qsh.cpp
  sto.f
  xem.f
  xpr.f
  xsc.f
  xse.f
  zzz.f
)

if(WITH_FASTJET)
  list(APPEND TP_SRC fastj.cpp)
else()
  list(APPEND TP_SRC x_j.f)
endif()

add_library(TP OBJECT ${TP_SRC})

if(WITH_FASTJET)
  target_include_directories(TP PRIVATE ${FASTJET_INCLUDE_DIR} ../KW)
else()
  target_include_directories(TP PRIVATE ../KW)
endif()
