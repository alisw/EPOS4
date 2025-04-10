#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# Do the local stuff:
set(KWi_SRC
  Accum.cpp
  Array.cpp
  Core.cpp
  Ccoti.cpp
  Eosp.cpp
  Event.cpp
  Files.cpp
  Hydyn.cpp
  Ico.cpp
  Outlist.cpp  
  Resc.cpp
)

add_library(KWi OBJECT ${KWi_SRC})
separate_arguments(CXX_COMPILE_FLAGS UNIX_COMMAND ${CXX_COMPILE_FLAGS})
target_compile_options(KWi PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CXX_COMPILE_FLAGS}>)


