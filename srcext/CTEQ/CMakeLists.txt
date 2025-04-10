# Do the local stuff:
if(WITH_CTEQ)
  set(CTEQ_SRC 
    CT14Pdf.f
    interface.f
    )
else()
  # HQ is unused... link with a dummy interface
  set(CTEQ_SRC
    dummyInterface.f
  )
  # Fortran compiler
  # Define this directory's flags:  
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64|aarch64")        
      list(APPEND FORTRAN_COMPILE_FLAGS "-Wno-unused-dummy-argument") 
      list(JOIN FORTRAN_COMPILE_FLAGS " " FORTRAN_COMPILE_FLAGS)
    else()
      message(ERROR "EPOS: EPOS compilation only succeeds on x86_64 or aarch64 architecture")
    endif()        
  else()
    message(ERROR "EPOS: Compilation with any compiler other than GCC has not yet been tested")
  endif()
endif()

add_library(CTEQ OBJECT ${CTEQ_SRC})
separate_arguments(FORTRAN_COMPILE_FLAGS UNIX_COMMAND ${FORTRAN_COMPILE_FLAGS})
target_compile_options(CTEQ PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAGS}>)
