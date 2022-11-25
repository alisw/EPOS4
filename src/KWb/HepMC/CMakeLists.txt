#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# Define this directory's flags:
set(${PROJECT_NAME}_HEPMC_CXX_FLAGS "-O2 -g -Wall -fPIC -pthread -m64" CACHE STRING "${PROJECT_NAME} HEPMC CXX FLAGS")
# Enable this directory's flags:
set(CMAKE_CXX_FLAGS "${${PROJECT_NAME}_HEPMC_CXX_FLAGS}")

# Do the local stuff:
set(HepMC_SRC 
  CompareGenEvent.cpp 
  PdfInfo.cpp 
  Flow.cpp 
  GenCrossSection.cpp 
  GenEvent.cpp 
  GenEventStreamIO.cpp 
  GenParticle.cpp
  GenRanges.cpp 
  GenVertex.cpp 
  HeavyIon.cpp 
  IO_AsciiParticles.cpp 
  IO_GenEvent.cpp 
  PdfInfo.cpp
  Polarization.cpp 
  SearchVector.cpp
  StreamHelpers.cpp
  StreamInfo.cpp
  Units.cpp
  WeightContainer.cpp 
  filterEvent.cpp 
)

# add_library(HepMC SHARED ${HepMC_SRC})
## set_target_properties(bg PROPERTIES LINKER_LANGUAGE FORTRAN)
## target_include_directories(bg PRIVATE ../KW/)
# install(TARGETS HepMC DESTINATION lib)

add_library(HEPMC OBJECT ${HepMC_SRC})
