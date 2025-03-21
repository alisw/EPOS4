#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# Do the local stuff:
set(YK_SRC 
  cll.cpp 
  conv.cpp  
  cornelius.cpp 
  ctrl.cpp
  dbs.cpp
  dkyc.cpp 
  es.cpp
  fld.cpp 
  hdo.cpp 
  ickw.cpp
  ihi.cpp 
  pdg.cpp 
  ptl.cpp 
  rmn.cpp  
  trancoeff.cpp
  uku.cpp 
)

add_library(YK OBJECT ${YK_SRC})
separate_arguments(CXX_COMPILE_FLAGS UNIX_COMMAND ${CXX_COMPILE_FLAGS})
target_compile_options(YK PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CXX_COMPILE_FLAGS}>)

target_include_directories(YK PRIVATE
  .
  ${PROJECT_SOURCE_DIR}/src/MS/
  ${ROOT_INCLUDE_DIRS}
)
