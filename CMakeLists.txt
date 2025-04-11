#
#  This file is part of EPOS4
#  Copyright (C) 2022 research institutions and authors (See CREDITS file)
#  This file is distributed under the terms of the GNU General Public License version 3 or later
#  (See COPYING file for the text of the licence)
#

# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 2.8.12.2...3.17.5)

# This is your project statement. You should always list languages.
project(EPOS LANGUAGES Fortran C CXX VERSION 4.0.0)

include("GNUInstallDirs")

option(INTERNAL_USE "internal use" FALSE)

# @ONLY : Restrict variable replacement to references of the form @VAR@. This is useful for configuring scripts that use ${VAR} syntax.
configure_file(scripts/epos.in scripts/epos @ONLY)

if(INTERNAL_USE)
  message(STATUS "EPOS: compilation in mode INTERNAL_USE")
  configure_file(scripts/epos+.in scripts/epos+ @ONLY)
  configure_file(scripts/epos+c.in scripts/epos+c @ONLY)
endif()

set(CMAKE_INSTALL_MESSAGE LAZY)
  
# Compilation option FULL, DEMO, BASIC or PHSD
if(NOT DEFINED COMPILE_OPTION)
  if(DEFINED ENV{COP})
    set(COMPILE_OPTION $ENV{COP})
  else()
    message(FATAL_ERROR "EPOS: COMPILE_OPTION is not defined")
  endif()
endif()

if(COMPILE_OPTION STREQUAL "FULL")
  set(WITH_ROOT True)
  set(WITH_BS False) 
  set(WITH_FASTJET True)
  set(WITH_HepMC3 True)
  set(WITH_HYDRODYNAMIC True)
  set(WITH_EB False)
  set(WITH_HQ True)
  set(WITH_CTEQ True)
  set(WITH_XA True)
  set(WITH_JT True)
elseif(COMPILE_OPTION STREQUAL "DEMO") 
  set(WITH_ROOT False)
  set(WITH_BS True) 
  set(WITH_FASTJET False)
  set(WITH_HepMC3 False)
  set(WITH_HYDRODYNAMIC False)
  set(WITH_EB False)
  set(WITH_HQ False)
  set(WITH_CTEQ False)
  set(WITH_XA False)
  set(WITH_JT False)
elseif(COMPILE_OPTION STREQUAL "BASIC") 
  set(WITH_ROOT True) 
  set(WITH_BS False) 
  set(WITH_FASTJET True)
  set(WITH_HepMC3 True)  
  set(WITH_HYDRODYNAMIC True)
  set(WITH_EB False)
  set(WITH_HQ False)
  set(WITH_CTEQ False)
  set(WITH_XA True)
  set(WITH_JT False)
elseif(COMPILE_OPTION STREQUAL "HQ") 
  set(WITH_ROOT True) 
  set(WITH_BS False) 
  set(WITH_FASTJET True)
  set(WITH_HepMC3 True)  
  set(WITH_HYDRODYNAMIC True)
  set(WITH_EB False)
  set(WITH_HQ True)
  set(WITH_CTEQ False)
  set(WITH_XA True)
  set(WITH_JT False)
elseif(COMPILE_OPTION STREQUAL "PHSD") 
  set(WITH_ROOT True) 
  set(WITH_BS True) 
  set(WITH_FASTJET False)
  set(WITH_HepMC3 False)
  set(WITH_HYDRODYNAMIC False)
  set(WITH_EB True)
  set(WITH_HQ False)
  set(WITH_CTEQ False)
  set(WITH_XA False)
  set(WITH_JT False)
endif()

message(STATUS "EPOS: Compilation with ${COMPILE_OPTION} option")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules" CACHE PATH "Module Path" FORCE)


if (NOT EPOS_CXX_STANDARD)
  set(EPOS_CXX_STANDARD 11)
endif()


if(WITH_ROOT)
  message(STATUS "EPOS: Looking for ROOT")
  find_package(ROOT 6.16 CONFIG REQUIRED COMPONENTS Tree Hist Physics)
  message(STATUS "EPOS: ROOT version ${ROOT_VERSION}")
  message(STATUS "EPOS: ROOT libraries ${ROOT_LIBRARIES}")
  message(STATUS "EPOS: ROOT include directories ${ROOT_INCLUDE_DIRS}")
  message(STATUS "EPOS: ROOT CXXFLAGS ${ROOT_CXX_FLAGS}")

  if(DEFINED ROOT_CXX_STANDARD)
    set(EPOS_CXX_STANDARD ${ROOT_CXX_STANDARD})
  else()
    string(REPLACE " " ";" ROOT_CXX_FLAGS_LIST ${ROOT_CXX_FLAGS})
    foreach(fl ${ROOT_CXX_FLAGS_LIST})
      if ("${fl}" MATCHES   ".*\\+\\+14$" OR "${fl}" MATCHES   ".*\\+\\+1y$" )
        set(EPOS_CXX_STANDARD 14)
      endif()
      if ("${fl}" MATCHES   ".*\\+\\+17$" OR "${fl}" MATCHES   ".*\\+\\+1z$" )
        set(EPOS_CXX_STANDARD 17)
      endif()
      if ("${fl}" MATCHES   ".*\\+\\+20$" OR "${fl}" MATCHES   ".*\\+\\+2a$" )
        set(EPOS_CXX_STANDARD 20)
      endif()
    endforeach(fl ${ROOT_CXX_FLAGS_LIST})
  endif()
endif()

message(STATUS "EPOS: EPOS_CXX_STANDARD : ${EPOS_CXX_STANDARD}")
set(CMAKE_CXX_STANDARD ${EPOS_CXX_STANDARD})

# FASTJET
if(WITH_FASTJET)
  message(STATUS "EPOS: Looking for FastJet")
  find_package(FastJet REQUIRED)
  message(STATUS "EPOS: FASTJET include directory ${FASTJET_INCLUDE_DIR}")
  message(STATUS "EPOS: FASTJET library ${FASTJET_LIBRARY}")
endif()

if(WITH_HepMC3)
  # ZLIB
  find_package(ZLIB REQUIRED)
  message(STATUS "EPOS: zlib found. zlib compressed files support enabled. ")
  message(STATUS "EPOS: ZLIB_INCLUDE_DIR : ${ZLIB_INCLUDE_DIR}")
  message(STATUS "EPOS: ZLIB_LIBRARIES : ${ZLIB_LIBRARIES}")
  
  # HEPMC3
  find_package(HepMC3 3.2.0 REQUIRED)
  message(STATUS "EPOS: HEPMC3_INCLUDE_DIR : ${HEPMC3_INCLUDE_DIR}")
  message(STATUS "EPOS: HEPMC3_LIB : ${HEPMC3_LIB}")
  message(STATUS "EPOS: HEPMC3_SEARCH_LIB : ${HEPMC3_SEARCH_LIB}")
  message(STATUS "EPOS: HEPMC3 version : ${HEPMC3_VERSION}")      
  message(STATUS "EPOS: HEPMC3 libraries : ${HEPMC3_LIBRARIES}")    
endif()

# Fortran module files are created in $LIBDIR/Mod                                          
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/FortranModules)


message(STATUS "EPOS: compilation options")
message(STATUS "EPOS: system processor: " ${CMAKE_SYSTEM_PROCESSOR})
message(STATUS "EPOS: Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
message(STATUS "EPOS: CXX: " ${CMAKE_CXX_COMPILER_ID})
message(STATUS "EPOS: C: " ${CMAKE_C_COMPILER_ID})
# Define the project's compilation flags:
# Fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64|aarch64")
    list(APPEND FORTRAN_COMPILE_FLAGS "-mcmodel=large")
    list(APPEND FORTRAN_COMPILE_FLAGS "-cpp")
    list(APPEND FORTRAN_COMPILE_FLAGS "-fbacktrace")
    list(APPEND FORTRAN_COMPILE_FLAGS "-ffpe-trap=invalid,zero")
    list(APPEND FORTRAN_COMPILE_FLAGS "-fno-automatic")
    list(APPEND FORTRAN_COMPILE_FLAGS "-fcheck=bounds")
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wall")
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wtabs")
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wno-conversion")
    list(APPEND FORTRAN_COMPILE_FLAGS "-Wuninitialized")
    list(APPEND FORTRAN_COMPILE_FLAGS "-fno-var-tracking-assignments")
    list(APPEND FORTRAN_COMPILE_FLAGS "-std=legacy")
    list(JOIN FORTRAN_COMPILE_FLAGS " " FORTRAN_COMPILE_FLAGS)
    message(STATUS "EPOS: FORTRAN_COMPILE_FLAGS: " ${FORTRAN_COMPILE_FLAGS})
  else()
    message(ERROR "EPOS: EPOS compilation only succeeds on x86_64 or aarch64 architecture")
  endif()    
else()
  message(ERROR "EPOS: Compilation with any compiler other than GCC has not yet been tested")
endif()

# CXX compiler
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64|aarch64")
    list(APPEND CXX_COMPILE_FLAGS "-Wall")
    list(APPEND CXX_COMPILE_FLAGS "-mcmodel=large")
    list(JOIN CXX_COMPILE_FLAGS " " CXX_COMPILE_FLAGS)
    message(STATUS "EPOS: CXX_COMPILE_FLAGS: " ${CXX_COMPILE_FLAGS})
  else()
    message(ERROR "EPOS: EPOS compilation only succeeds on x86_64 or aarch64 architecture")
  endif()
else()
  message(ERROR "EPOS: Compilation with any compiler other than GCC has not yet been tested")
endif()

# C compiler
if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64|aarch64")
    list(APPEND C_COMPILE_FLAGS "-Wall")
    list(APPEND C_COMPILE_FLAGS "-mcmodel=large")
    list(JOIN C_COMPILE_FLAGS " " C_COMPILE_FLAGS)
    message(STATUS "EPOS: C_COMPILE_FLAGS: " ${C_COMPILE_FLAGS})
  else()
    message(ERROR "EPOS: EPOS compilation only succeeds on x86_64 or aarch64 architecture")    
  endif()
else()
  message(ERROR "EPOS: Compilation with any compiler other than GCC has not yet been tested")
endif()

## Enter the subdirectories from src:
add_subdirectory(interfaces)
add_subdirectory(src/KW)
add_subdirectory(src/KWa)
add_subdirectory(src/KWb)
add_subdirectory(src/KWi)
add_subdirectory(src/TP)
add_subdirectory(src/XA)

## Enter the subdirectories from srcext:
add_subdirectory(srcext/CTEQ) 
add_subdirectory(srcext/HQ)
add_subdirectory(srcext/XA)

if(WITH_EB)
  add_subdirectory(srcext/EB/phsd)
  add_subdirectory(srcext/EB)
else()
  add_subdirectory(srcext/UR)
endif()

if(WITH_JT)
  add_subdirectory(srcext/JT)
endif()

if(WITH_HYDRODYNAMIC)
  add_subdirectory(src/MS)
  add_subdirectory(srcext/YK)
endif()

####################################################################################
####################################################################################
# Xepos target
####################################################################################
####################################################################################
set(OBJECT_LIST
  $<TARGET_OBJECTS:INTERFACES> 
  $<TARGET_OBJECTS:KW>
  $<TARGET_OBJECTS:KWa> 
  $<TARGET_OBJECTS:KWb> 
  $<TARGET_OBJECTS:KWi> 
  $<TARGET_OBJECTS:TP> 
  $<TARGET_OBJECTS:CTEQ>	
  $<TARGET_OBJECTS:HQ>
  $<TARGET_OBJECTS:XA> 
  $<TARGET_OBJECTS:XAEXT> 
  )

if(WITH_JT)
  list(APPEND OBJECT_LIST
    $<TARGET_OBJECTS:JT> 
    )
endif()

if(WITH_HYDRODYNAMIC)
  list(APPEND OBJECT_LIST
    $<TARGET_OBJECTS:MS> 
    $<TARGET_OBJECTS:YK> 
    )
endif()

if(WITH_EB)
  list(APPEND OBJECT_LIST
    $<TARGET_OBJECTS:EB> 
    )
else()
  list(APPEND OBJECT_LIST
    $<TARGET_OBJECTS:UR> 
    )	
endif()

if(WITH_HQ)
  list(APPEND OBJECT_LIST
    $<TARGET_OBJECTS:HQ_Module> 
    )
endif()

if(WITH_XA)
  list(APPEND OBJECT_LIST
    $<TARGET_OBJECTS:XA> 
    )
endif()


message(STATUS "EPOS: OBJECT_LIST ${OBJECT_LIST}")

add_executable(Xepos src/main.cpp ${OBJECT_LIST})
set_target_properties(Xepos PROPERTIES LINK_FLAGS "-Wl,--no-relax")
set_target_properties(Xepos PROPERTIES ENABLE_EXPORTS ON)
set_target_properties(Xepos PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(WITH_HepMC3)
  target_link_libraries(Xepos PUBLIC ZLIB::ZLIB)  
  target_link_libraries(Xepos PUBLIC HepMC3::HepMC3)  
endif()
if(WITH_ROOT)
  target_link_libraries(Xepos PUBLIC ROOT::Tree ROOT::Hist ROOT::Physics)  
endif()
if(WITH_FASTJET)
  target_link_libraries(Xepos PUBLIC FastJet::FastJet)  
endif()
if(WITH_EB)
  target_link_libraries(Xepos PUBLIC PHSD EB_WRAPPER)
endif()

set(LINK_LIBRARY_LIST "")

list(APPEND LINK_LIBRARY_LIST
    $ENV{GMP_ROOT}/lib/libgmp.a
    )
list(APPEND LINK_LIBRARY_LIST
    $ENV{CGAL_ROOT}/lib/libCGAL.so
    )    

target_link_libraries(Xepos PUBLIC ${LINK_LIBRARY_LIST})


####################################################################################
####################################################################################
# tests
####################################################################################
####################################################################################
# testing binary
add_executable(Array_test test/Array_test.cpp src/KWi/Array.cpp)
target_include_directories(Array_test PRIVATE src/KWi)

# enable testing functionality
set(EPOS_ENABLE_TESTING ON)
if (EPOS_ENABLE_TESTING)
  ENABLE_TESTING()

  # see: https://klaus.pages.in2p3.fr/epos4/pdf/optnsExamples.pdf
  set(test_name test1 CACHE STRING "" FORCE)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name})
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name}/examples)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name}/histo)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name}/check)  
  file(COPY examples/proton_proton_scattering_no_hydro.optns
    DESTINATION  ${CMAKE_BINARY_DIR}/test1/examples/)
  add_test(NAME proton_proton_scattering_no_hydro
    COMMAND ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/epos proton_proton_scattering_no_hydro
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name}/)
  set_tests_properties(proton_proton_scattering_no_hydro PROPERTIES ENVIRONMENT 
    "EPO=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/;OPT=${CMAKE_BINARY_DIR}/${test_name}/examples/;HTO=${CMAKE_BINARY_DIR}/${test_name}/histo/;CHK=${CMAKE_BINARY_DIR}/${test_name}/check/")

  # see: https://klaus.pages.in2p3.fr/epos4/pdf/optnsExamples.pdf
  set(test_name test2 CACHE STRING "" FORCE)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name})
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name}/examples)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name}/histo)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name}/check)  
  file(COPY examples/proton_proton_scattering_with_hydro.optns
    DESTINATION  ${CMAKE_BINARY_DIR}/test2/examples/)
  add_test(NAME proton_proton_scattering_with_hydro
    COMMAND ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/epos proton_proton_scattering_with_hydro
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/${test_name}/)
  set_tests_properties(proton_proton_scattering_with_hydro PROPERTIES ENVIRONMENT 
    "EPO=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/;OPT=${CMAKE_BINARY_DIR}/${test_name}/examples/;HTO=${CMAKE_BINARY_DIR}/${test_name}/histo/;CHK=${CMAKE_BINARY_DIR}/${test_name}/check/")

# define tests
  add_test(
    NAME Array_test
    COMMAND $<TARGET_FILE:Array_test>
  )
endif()


####################################################################################
####################################################################################
# installation
####################################################################################
####################################################################################
install(TARGETS Xepos DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES ${CMAKE_BINARY_DIR}/scripts/epos
  DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT config
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
  GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

install(FILES ${CMAKE_SOURCE_DIR}/README.md ${CMAKE_SOURCE_DIR}/VERSION.txt
  DESTINATION ${CMAKE_INSTALL_DATADIR})
if(EXISTS ${CMAKE_SOURCE_DIR}/banner.txt)
install(FILES ${CMAKE_SOURCE_DIR}/banner.txt DESTINATION ${CMAKE_INSTALL_DATADIR})
endif()
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/KWt/ DESTINATION ${CMAKE_INSTALL_DATADIR}/src/KWt
  FILES_MATCHING PATTERN "*.dt" PATTERN "*.eos.gz" PATTERN "*.eos")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/KWn/ DESTINATION ${CMAKE_INSTALL_DATADIR}/src/KWn
  FILES_MATCHING PATTERN "*.i")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/srcext/UR/ DESTINATION ${CMAKE_INSTALL_DATADIR}/srcext/UR
  FILES_MATCHING PATTERN "*.dat")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/srcext/URt/ DESTINATION ${CMAKE_INSTALL_DATADIR}/srcext/URt
  FILES_MATCHING PATTERN "*.dat")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/srcext/YK/ DESTINATION ${CMAKE_INSTALL_DATADIR}/srcext/YK
  FILES_MATCHING PATTERN "*.data")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/config/ DESTINATION ${CMAKE_INSTALL_DATADIR}/config
  FILES_MATCHING PATTERN "*.i")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ DESTINATION ${CMAKE_INSTALL_DATADIR}/examples
  FILES_MATCHING PATTERN "*.optns")
if(WITH_HQ)
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/srcext/HQt/ DESTINATION ${CMAKE_INSTALL_DATADIR}/srcext/HQt
    FILES_MATCHING PATTERN "*.dat")
endif()  

####################################################################################
####################################################################################
# XeposC target
####################################################################################
####################################################################################
if(INTERNAL_USE)
  add_executable(XeposC $<TARGET_OBJECTS:XEPOSC>)
  set_target_properties(XeposC PROPERTIES LINKER_LANGUAGE CXX)
  set_target_properties(XeposC PROPERTIES ENABLE_EXPORTS ON)
  set_target_properties(XeposC PROPERTIES POSITION_INDEPENDENT_CODE ON)
  
  install(TARGETS XeposC DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${CMAKE_BINARY_DIR}/scripts/epos+
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT config
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
  install(FILES ${CMAKE_BINARY_DIR}/scripts/epos+c
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT config
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
endif()
