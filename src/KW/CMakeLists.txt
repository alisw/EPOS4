#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# Do the local stuff:
set(KW_SRC
  bas.f
  bas.cpp
  bor.f
  conf.f
  d.f
  dky.f
  ems.f
  fra.f
  freeze.f
  graph2.f
  graph.f
  h.f
  hg.f
  hnb.f
  ho.h
  ico.f
  ids.f
  int.f
  jpo.f
  key1.f90
  key2.f
  lea.f
  ${PROJECT_SOURCE_DIR}/src/modif.cpp
  ${PROJECT_SOURCE_DIR}/src/list.cpp  
  omg.f
  rsh.f
  sem.f
  tab.f
  tim.f
  uc.f
  uti.f
  xan.f
  xht.f
)

if(WITH_ROOT)
  list(APPEND KW_SRC
    out.f)
endif()

if(NOT WITH_EB)
  list(APPEND KW_SRC
    hyt.f
    eos.f
    u.f)
endif()

# Enable this directory's FORTRAN flags:
if(WITH_ROOT)
  list(APPEND FORTRAN_COMPILE_FLAGS "-D__ROOT__")
  list(JOIN FORTRAN_COMPILE_FLAGS " " FORTRAN_COMPILE_FLAGS)
endif()

if(WITH_BS)
  list(APPEND FORTRAN_COMPILE_FLAGS "-D__BS__")
  list(JOIN FORTRAN_COMPILE_FLAGS " " FORTRAN_COMPILE_FLAGS)
endif()

if(WITH_EB)
  list(APPEND FORTRAN_COMPILE_FLAGS "-D__EB__")
  list(JOIN FORTRAN_COMPILE_FLAGS " " FORTRAN_COMPILE_FLAGS)
endif()


add_library(KW OBJECT ${KW_SRC})
separate_arguments(FORTRAN_COMPILE_FLAGS UNIX_COMMAND ${FORTRAN_COMPILE_FLAGS})
target_compile_options(KW PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAGS}>)
target_include_directories(
    KW
    PRIVATE ${PROJECT_SOURCE_DIR}/src/TP/
    PRIVATE ${PROJECT_SOURCE_DIR}/src/KWa/
    PRIVATE ${PROJECT_SOURCE_DIR}/src/KWb/
    PRIVATE ${PROJECT_SOURCE_DIR}/src/KWi/
    PRIVATE ${PROJECT_SOURCE_DIR}/srcext/UR/)

if(WITH_ROOT)
  target_link_libraries(KW PRIVATE ROOT::Core)  
endif()


if(INTERNAL_USE)
  # Target XeposC
  set(XEPOSC_SRC
    c.f
    d.f
  )

  add_library(XEPOSC OBJECT ${XEPOSC_SRC})
  target_compile_options(XEPOSC PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAGS}>)
endif()
