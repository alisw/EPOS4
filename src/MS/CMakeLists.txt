#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# Do the local stuff:
set(MS_SRC 
  eo1.cpp
  eo3.cpp
  eoBEST.cpp   
  eoCEM.cpp
  eoChiral.cpp
  eoPNJL.cpp   
  eos.cpp
)

add_library(MS OBJECT ${MS_SRC})
separate_arguments(CXX_COMPILE_FLAGS UNIX_COMMAND ${CXX_COMPILE_FLAGS})
target_compile_options(MS PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CXX_COMPILE_FLAGS}>)

target_include_directories(MS PRIVATE ${PROJECT_SOURCE_DIR}/srcext/YK/)
