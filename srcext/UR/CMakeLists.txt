#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#


# Define the compilation flags:
# Fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64|aarch64")      
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wno-unused-variable")
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wno-uninitialized")
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wno-unused-dummy-argument")
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wno-unused-label")
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wno-integer-division")
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wno-character-truncation")
    list(JOIN FORTRAN_COMPILE_FLAGS " " FORTRAN_COMPILE_FLAGS)
  else()
    message(ERROR "EPOS: EPOS compilation only succeeds on x86_64 or aarch64 architecture")
  endif()        
else()
  message(ERROR "EPOS: Compilation with any compiler other than GCC has not yet been tested")
endif()

# Do the local stuff:
set(UR_SRC 
  addpart.f  
  angdis.f   
  anndec.f 
  blockres.f 
# boxinc.f 
  boxprg.f   
  cascinit.f   
# colltab.f   
  coload.f   
# comnorm.f 
# comres.f 
# coms.f      
# comstr.f  
# comwid.f 
  dectim.f   
  delpart.f    
  detbal.f   
  dwidth.f   
  error.f      
# freezeout.f 
  getmass.f  
  getspin.f  
  hepchg.f     
  hepcmp.f   
  hepnam.f   
  init.f       
  input.f 
# inputs.f  
  iso.f      
  ityp2pdg.f 
  jdecay2.f    
  make22.f   
# newpart.f 
  numrec.f   
# options.f  
# outcom.f
  output.f     
  paulibl.f 
  proppot.f  
  saveinfo.f 
  scatter.f    
  siglookup.f 
  string.f  
  tabinit.f    
  uaaa.f     
  upmerge.f  
  upyth.f
  whichres.f   
)

add_library(UR OBJECT ${UR_SRC})
separate_arguments(FORTRAN_COMPILE_FLAGS UNIX_COMMAND ${FORTRAN_COMPILE_FLAGS})
target_compile_options(UR PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAGS}>)
