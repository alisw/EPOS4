#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#


if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64|aarch64")  
    list(APPEND CXX_COMPILE_FLAGS "-Wno-sign-compare")
    list(APPEND CXX_COMPILE_FLAGS "-Wno-unused-variable")
    list(APPEND CXX_COMPILE_FLAGS "-Wno-maybe-uninitialized")
    list(APPEND CXX_COMPILE_FLAGS "-Wno-reorder")
    list(JOIN CXX_COMPILE_FLAGS " " CXX_COMPILE_FLAGS)
  else()
    message(ERROR "EPOS: EPOS compilation only succeeds on x86_64 or aarch64 architecture")
  endif()        
else()
  message(ERROR "EPOS: Compilation with any compiler other than GCC has not yet been tested")
endif()

# Do the local stuff:
set(KWa_SRC
  Subanalysis.cpp  
  Analysis.cpp
  AnalysisInterface.cpp
  Utils.cpp
)

add_library(KWa OBJECT ${KWa_SRC})
separate_arguments(CXX_COMPILE_FLAGS UNIX_COMMAND ${CXX_COMPILE_FLAGS})
target_compile_options(KWa PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CXX_COMPILE_FLAGS}>)

target_include_directories(KWa PRIVATE
  ${PROJECT_SOURCE_DIR}/src/KWa
  ${PROJECT_SOURCE_DIR}/srcext/XA/include
)
