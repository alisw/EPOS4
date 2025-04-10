#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#
# Define this directory's flags:
if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64|aarch64")    
    list(APPEND C_COMPILE_FLAGS "-Wno-deprecated-declarations")
    list(APPEND C_COMPILE_FLAGS "-Wno-unused-variable")
    list(APPEND C_COMPILE_FLAGS "-Wno-unused-result")
    list(JOIN C_COMPILE_FLAGS " " C_COMPILE_FLAGS)
  else()
    message(ERROR "EPOS: EPOS compilation only succeeds on x86_64 or aarch64 architecture")
  endif()        
else()
    message(ERROR "EPOS: Compilation with any compiler other than GCC has not yet been tested")
endif()
list(APPEND CXX_COMPILE_FLAGS "-DHEPMC3_Z_SUPPORT")
list(APPEND CXX_COMPILE_FLAGS "-DHEPMC3_USE_COMPRESSION")
list(JOIN CXX_COMPILE_FLAGS " " CXX_COMPILE_FLAGS)

if(WITH_ROOT)
  list(APPEND CXX_COMPILE_FLAGS "-D__ROOT__")
  list(JOIN CXX_COMPILE_FLAGS " " CXX_COMPILE_FLAGS)  
endif()

# Do the local stuff:
if(WITH_ROOT)
  set(INTERFACES_SRC
    EposRoot.cpp
    EposHepmc.cpp
  )
endif()

if(WITH_ROOT)
  add_library(INTERFACES OBJECT ${INTERFACES_SRC})
  separate_arguments(CXX_COMPILE_FLAGS UNIX_COMMAND ${CXX_COMPILE_FLAGS})
  target_compile_options(INTERFACES PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CXX_COMPILE_FLAGS}>)     
  target_include_directories(INTERFACES PRIVATE ${ROOT_INCLUDE_DIRS} ${HEPMC3_INCLUDE_DIR})
endif()
