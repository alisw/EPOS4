# find_package(ROOT)
# Define this directory's flags:
# set(${PROJECT_NAME}_YK_CXX_FLAGS "-std=c++11 -O2 -g -Wall -fPIC" CACHE STRING "${PROJECT_NAME} YK CXX FLAGS")

list(APPEND CXX_FLAGS ${ROOT_CPPFLAGS})
list(APPEND CXX_FLAGS "-O2 -g -Wall -fPIC")
list(JOIN CXX_FLAGS " " ${PROJECT_NAME}_YK_CXX_FLAGS)


# Enable this directory's flags:
set(CMAKE_CXX_FLAGS "${${PROJECT_NAME}_YK_CXX_FLAGS}")
message(STATUS "YK CXX FLAGS: " ${CMAKE_CXX_FLAGS})   

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
target_include_directories(YK PRIVATE ../../src/MS/ PRIVATE . PRIVATE ${ROOT_INCLUDE_DIR})
