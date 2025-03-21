# - Locate FastJet library
# Defines:
#
#  FASTJET_FOUND
#  FASTJET_INCLUDE_DIR
#  FASTJET_INCLUDE_DIRS (not cached)
#  FASTJET_LIBRARY
#  FASTJET_LIBRARIES (not cached)
#  FASTJET_LIBRARY_DIRS (not cached)

find_path(FASTJET_INCLUDE_DIR fastjet/version.hh
          HINTS /usr/include /usr/local/include $ENV{FASTSYS}/include ${FASTSYS}/include)

find_library(FASTJET_LIBRARY NAMES fastjet
             HINTS /usr/local/lib /usr/lib /usr/lib64 $ENV{FASTSYS}/lib ${FASTSYS}/lib)

# handle the QUIETLY and REQUIRED arguments and set FASTJET_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FastJet DEFAULT_MSG FASTJET_INCLUDE_DIR FASTJET_LIBRARY)

mark_as_advanced(FASTJET_FOUND FASTJET_INCLUDE_DIR FASTJET_LIBRARY)

if(FASTJET_FOUND AND NOT TARGET FastJet::FastJet)
    add_library(FastJet::FastJet UNKNOWN IMPORTED)
    set_target_properties(FastJet::FastJet PROPERTIES
        IMPORTED_LOCATION "${FASTJET_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FASTJET_INCLUDE_DIR}"
    )
endif()

