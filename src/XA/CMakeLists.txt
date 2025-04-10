#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

if(WITH_XA)
  set(XA_SRC 
    z1.f
    z2.f 
    z4.c
    )
else()
  # XA is unused... link with a dummy interface
  set(XA_SRC
    dummyInterface.f
    dummyInterface.c
    )
    # Define this directory's flags:
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
      if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64|aarch64")          
	list(APPEND FORTRAN_COMPILE_FLAGS "-Wno-unused-dummy-argument") 
	list(JOIN FORTRAN_COMPILE_FLAGS " " FORTRAN_COMPILE_FLAGS)
      else()
	message(ERROR "EPOS: EPOS compilation only succeeds on x86_64 or aarch64 architecture")
      endif()        
    else()
      message(ERROR "EPOS: Compilation with any compiler other than GCC has not yet been tested")
    endif()
endif()

add_library(XA OBJECT ${XA_SRC})
separate_arguments(C_COMPILE_FLAGS UNIX_COMMAND ${C_COMPILE_FLAGS})
separate_arguments(FORTRAN_COMPILE_FLAGS UNIX_COMMAND ${FORTRAN_COMPILE_FLAGS})
target_compile_options(XA PUBLIC $<$<COMPILE_LANGUAGE:C>:${C_COMPILE_FLAGS}>)
target_compile_options(XA PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAGS}>)

#if(WITH_XA)
  # add the module dependency if nedded
#  add_dependencies(XA HQ_Module)
#endif()

target_include_directories(XA PRIVATE
  ${PROJECT_SOURCE_DIR}/src/KW/
  ${PROJECT_SOURCE_DIR}/src/KWa/
)
