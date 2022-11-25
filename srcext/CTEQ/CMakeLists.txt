# Do the local stuff:
if(WITH_CTEQ)
  set(CTEQ_SRC 
    CT14Pdf.f
    interface.f
    )
  # Define this directory's flags:
  set(${PROJECT_NAME}_CTEQ_FORTRAN_FLAGS "-g -mcmodel=large -cpp -fbacktrace -ffpe-trap=invalid,zero  -Wall -fno-automatic -fcheck=bounds -O2 -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments" CACHE STRING "${PROJECT_NAME} CTEQ FORTRAN FLAGS")
else()
  # VO is unused... link with a dummy interface
  set(CTEQ_SRC
    dummyInterface.f
    )
  # Define this directory's flags:
  set(${PROJECT_NAME}_CTEQ_FORTRAN_FLAGS "-g -mcmodel=large -fPIC -cpp -fbacktrace -ffpe-trap=invalid,zero -fno-automatic -fcheck=bounds -O2 -Wall -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments -Wno-unused-dummy-argument" CACHE STRING "${PROJECT_NAME} CTEQ FORTRAN FLAGS")
endif()

# Enable this directory's flags:
set(CMAKE_Fortran_FLAGS "${${PROJECT_NAME}_CTEQ_FORTRAN_FLAGS}")

#add_library(CTEQ SHARED ${CTEQ_SRC})
## set_target_properties(bg PROPERTIES LINKER_LANGUAGE FORTRAN)
## target_include_directories(CTEQ PRIVATE ../KW/)
#install(TARGETS CTEQ DESTINATION lib)
add_library(CTEQ OBJECT ${CTEQ_SRC})
