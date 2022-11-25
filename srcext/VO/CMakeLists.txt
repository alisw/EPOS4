# Do the local stuff:
if(WITH_VO)
  # VO is used
  set(VO_SRC 
    main7.f 
    pack0.f 
    pack1.f 
    pack2.f 
    pack3.f 
    pack4.f 
    pack5.f 
    pack6.f 
    pack7.f 
    pack8.f 
    packboltzmann_radiat_v3.f 
    packfraganddecay.f 
    packspectra.f
    packeposhq.f
    packeposhq.cpp
    )
  
  # Define this directory's flags:
  set(${PROJECT_NAME}_VO_FORTRAN_FLAGS "-g -mcmodel=large -fPIC -cpp -fbacktrace -ffpe-trap=invalid,zero -fno-automatic -fcheck=bounds -O2 -Wall -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments" CACHE STRING "${PROJECT_NAME} VO FORTRAN FLAGS")

else()
  # VO is unused... link with a dummy interface
  set(VO_SRC
    dummyInterface.f
    packeposhq.cpp
    )
  # Define this directory's flags:
  set(${PROJECT_NAME}_VO_FORTRAN_FLAGS "-g -mcmodel=large -fPIC -cpp -fbacktrace -ffpe-trap=invalid,zero -fno-automatic -fcheck=bounds -O2 -Wall -Wtabs -Wno-conversion -Wuninitialized -fno-var-tracking-assignments -Wno-unused-dummy-argument" CACHE STRING "${PROJECT_NAME} VO FORTRAN FLAGS")
endif()

# Enable this directory's flags:
set(CMAKE_Fortran_FLAGS "${${PROJECT_NAME}_VO_FORTRAN_FLAGS}")

# Define this directory's CXX flags:
set(${PROJECT_NAME}_VO_CXX_FLAGS "-O2 -g -Wall -fPIC -pthread -m64" CACHE STRING "${PROJECT_NAME} VO CXX FLAGS")
# Enable this directory's CXX flags:
set(CMAKE_CXX_FLAGS "${${PROJECT_NAME}_VO_CXX_FLAGS}")

add_library(VO OBJECT ${VO_SRC})
target_include_directories(VO PRIVATE ../../src/KW/)
