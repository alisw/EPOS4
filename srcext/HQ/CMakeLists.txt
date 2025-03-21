#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# Do the local stuff:
if(WITH_HQ)
  # compile the fortran modules first
  set(HQ_Module_SRC
    main7.f
  )
  add_library(HQ_Module OBJECT ${HQ_Module_SRC})

  # HQ is used
  set(HQ_SRC 
    pack0.f 
    pack1.f 
    pack2.f 
    pack3.f 
    pack4.f 
    pack5.f 
    pack6.f 
    pack7.f 
    pack8.f
    pack9.f
    packboltzmann_radiat_v3.f 
    packfraganddecay.f 
    packspectra.f
    packeposhq.f
    packeposhq.cpp
    )
else()
  # HQ is unused... link with a dummy interface
  set(HQ_SRC
    dummyInterface.f
    packeposhq.cpp
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

add_library(HQ OBJECT ${HQ_SRC})
separate_arguments(FORTRAN_COMPILE_FLAGS UNIX_COMMAND ${FORTRAN_COMPILE_FLAGS})
target_compile_options(HQ PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAGS}>)

if(WITH_HQ)
  # add the module dependency if nedded
  add_dependencies(HQ HQ_Module)
endif()

target_include_directories(HQ PRIVATE ${PROJECT_SOURCE_DIR}/src/KW/)


