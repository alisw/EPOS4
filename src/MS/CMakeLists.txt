# Define this directory's flags:
set(${PROJECT_NAME}_MS_CXX_FLAGS "-O2 -g -Wall -fPIC -pthread -m64" CACHE STRING "${PROJECT_NAME} HEPMC CXX FLAGS")
# Enable this directory's flags:
set(CMAKE_CXX_FLAGS "${${PROJECT_NAME}_MS_CXX_FLAGS}")

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
target_include_directories(MS PRIVATE ../../srcext/YK/)
