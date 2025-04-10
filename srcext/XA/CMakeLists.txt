#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

set(XA_SRCEXT
  PtSpectra.cpp
  PtSpectraMultTrigger.cpp
  MyAnalysis.cpp
  ParticlePairsCentralityTrigger.cpp
  ParticleSpectraMultiplicityTrigger.cpp
  ParticleSpectraCentralityTrigger.cpp
  ParticleSpectraPercentileTrigger.cpp
  MultiplicityDistributions.cpp
  MultiplicityDependences.cpp
  FlowScalarProductCentralityTrigger.cpp
  )

add_library(XAEXT OBJECT ${XA_SRCEXT})
separate_arguments(CXX_COMPILE_FLAGS UNIX_COMMAND ${CXX_COMPILE_FLAGS})
target_compile_options(XAEXT PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CXX_COMPILE_FLAGS}>)

target_include_directories(XAEXT PRIVATE
  include/
  ${PROJECT_SOURCE_DIR}/src/KW/
  ${PROJECT_SOURCE_DIR}/src/KWa/
)
