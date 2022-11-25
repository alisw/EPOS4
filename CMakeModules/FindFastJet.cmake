#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

#####################################################################################
# (c) Copyright 1998-2019 CERN for the benefit of the LHCb and ATLAS collaborations #
#                                                                                   #
# This software is distributed under the terms of the Apache version 2 licence,     #
# copied verbatim in the file "LICENSE".                                            #
#                                                                                   #
# In applying this licence, CERN does not waive the privileges and immunities       #
# granted to it by virtue of its status as an Intergovernmental Organization        #
# or submit itself to any jurisdiction.                                             #
#####################################################################################
# - Locate FastJet library
# Defines:
#
#  FASTJET_FOUND
#  FASTJET_INCLUDE_DIR
#  FASTJET_INCLUDE_DIRS (not cached)
#  FASTJET_LIBRARY
#  FASTJET_LIBRARIES (not cached)
#  FASTJET_LIBRARY_DIRS (not cached)


# the following disables all default paths (either from cmake, from environment)
FIND_PATH (FASTJET_BIN_DIR fastjet-config
  PATHS 
  ${FASTJETSYS}/bin	
  $ENV{FASTJETSYS}/bin
  $ENV{AUGER_BASE}/External/FASTJET/pro/bin
  NO_DEFAULT_PATH
)

# now defaults are allowed but if nothing is found it is not overwritten
# (because it is cached)
FIND_PATH (FASTJET_BIN_DIR fastjet-config)

GET_FILENAME_COMPONENT (FASTJETSYS ${FASTJET_BIN_DIR} PATH)

IF (NOT ENV{FASTJETSYS})
  SET (ENV{FASTJETSYS} ${FASTJETSYS})
ENDIF (NOT ENV{FASTJETSYS})

IF (FASTJETSYS)
  # ----------------------------------------------------------------------------
  # Get ROOT version, re-express in form XX.YY.ZZ, compare to requirements.
  # ----------------------------------------------------------------------------
  EXECUTE_PROCESS (COMMAND ${FASTJET_BIN_DIR}/fastjet-config --version
    OUTPUT_VARIABLE FASTJET_VERSION
    )
  STRING (REGEX REPLACE "[ \t\r\n]+" "" FASTJET_VERSION "${FASTJET_VERSION}") 
  
  find_path(FASTJET_INCLUDE_DIR fastjet/version.hh
    HINTS $ENV{FASTJET_ROOT_DIR}/include ${FASTJET_ROOT_DIR}/include)
  
  find_library(FASTJET_LIBRARY NAMES fastjet
    HINTS $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)
  
  # handle the QUIETLY and REQUIRED arguments and set FASTJET_FOUND to TRUE if
  # all listed variables are TRUE
  INCLUDE(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(FastJet DEFAULT_MSG FASTJET_INCLUDE_DIR FASTJET_LIBRARY)
  
  mark_as_advanced(FASTJET_FOUND FASTJET_INCLUDE_DIR FASTJET_LIBRARY)
  
  set(FASTJET_INCLUDE_DIRS ${FASTJET_INCLUDE_DIR})
  set(FASTJET_LIBRARIES ${FASTJET_LIBRARY})
  get_filename_component(FASTJET_LIBRARY_DIRS ${FASTJET_LIBRARY} PATH)

ENDIF (FASTJETSYS)
