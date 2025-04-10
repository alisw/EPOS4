#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64|aarch64")      
    list(APPEND CXX_COMPILE_FLAGS "-Wno-deprecated-declarations")
    list(JOIN CXX_COMPILE_FLAGS " " CXX_COMPILE_FLAGS)
  else()
    message(ERROR "EPOS: EPOS compilation only succeeds on x86_64 or aarch64 architecture")
  endif()            
else()
  message(ERROR "EPOS: Compilation with any compiler other than GCC has not yet been tested")
endif()
  

# Do the local stuff:
set(TP_SRC 
  gen.f
  models.f  
  par.f 
  pdf.f
  ptb.f
  qsh.f
  qsh.cpp
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
separate_arguments(CXX_COMPILE_FLAGS UNIX_COMMAND ${CXX_COMPILE_FLAGS})
separate_arguments(FORTRAN_COMPILE_FLAGS UNIX_COMMAND ${FORTRAN_COMPILE_FLAGS})
target_compile_options(TP PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CXX_COMPILE_FLAGS}>)
target_compile_options(TP PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAGS}>)

if(WITH_FASTJET)
  target_include_directories(TP PRIVATE
    ${FASTJET_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/src/KW
  )
else()
  target_include_directories(TP PRIVATE
    ${PROJECT_SOURCE_DIR}/src/KW)
endif()
