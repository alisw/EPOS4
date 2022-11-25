# Define this directory's flags:
set(${PROJECT_NAME}_UR_FORTRAN_FLAGS "-g -mcmodel=large -fPIC -cpp -fbacktrace -ffpe-trap=invalid,zero -fno-automatic -fcheck=bounds" CACHE STRING "${PROJECT_NAME} UR FORTRAN FLAGS")
# Enable this directory's flags:
set(CMAKE_Fortran_FLAGS "${${PROJECT_NAME}_UR_FORTRAN_FLAGS}")

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

#add_library(UR SHARED ${UR_SRC})
## set_target_properties(bg PROPERTIES LINKER_LANGUAGE FORTRAN)
## target_include_directories(UR PRIVATE ../KW/)
#install(TARGETS UR DESTINATION lib)

add_library(UR OBJECT ${UR_SRC})
