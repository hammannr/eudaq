if(NOT EUDAQ_BUILD_EXECUTABLE)
  message(STATUS "Disable the building of main EUDAQ executables (EUDAQ_BUILD_EXECUTABLE=OFF)")
  return()
endif()

set(EXE_CLI_TLU_READER euCliTluReader)
add_executable(${EXE_CLI_TLU_READER} src/euCliTluReader.cxx)
target_link_libraries(${EXE_CLI_TLU_READER} ${EUDAQ_CORE_LIBRARY} ${EUDAQ_THREADS_LIB})
list(APPEND INSTALL_TARGETS ${EXE_CLI_TLU_READER})

set(EXE_CLI_TRIGGER_READER euCliTriggerReader)
add_executable(${EXE_CLI_TRIGGER_READER} src/euCliTriggerReader.cxx)
target_link_libraries(${EXE_CLI_TRIGGER_READER} ${EUDAQ_CORE_LIBRARY} ${EUDAQ_THREADS_LIB})
list(APPEND INSTALL_TARGETS ${EXE_CLI_TRIGGER_READER})

if(USER_TLU_BUILD_EUDET)
  message(STATUS "Building EUDET TLU stand-alone executables (USER_BUILD_EUDET_TLU=ON)")

  set(EXE_EUDET_TLU_CONTROL EudetTluControl)
  add_executable(${EXE_EUDET_TLU_CONTROL} src/EudetTluControl.cxx)
  target_link_libraries(${EXE_EUDET_TLU_CONTROL} ${EUDAQ_CORE_LIBRARY} ${USER_HARDWARE_LIBRARY})
  list(APPEND INSTALL_TARGETS ${EXE_EUDET_TLU_CONTROL})

  set(EXE_EUDET_TLU_NOROOT EudetTluNoRoot)
  add_executable(${EXE_EUDET_TLU_NOROOT} src/EudetTluNoRoot.cxx)
  target_link_libraries(${EXE_EUDET_TLU_NOROOT} ${EUDAQ_CORE_LIBRARY} ${USER_HARDWARE_LIBRARY})
  list(APPEND INSTALL_TARGETS ${EXE_EUDET_TLU_NOROOT})

  set(EXE_EUDET_TLU_RESET EudetTluReset)
  add_executable(${EXE_EUDET_TLU_RESET} src/EudetTluReset.cxx)
  target_link_libraries(${EXE_EUDET_TLU_RESET} ${EUDAQ_CORE_LIBRARY} ${USER_HARDWARE_LIBRARY})
  list(APPEND INSTALL_TARGETS ${EXE_EUDET_TLU_RESET})
else()
  message(STATUS "EUDET TLU stand-alone executables are NOT to be built (USER_BUILD_EUDET_TLU=OFF)")
endif()

if(USER_TLU_BUILD_AIDA)
  message(STATUS "Building AIDA TLU stand-alone executables (USER_BUILD_AIDA_TLU=ON)")
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
  find_package(ROOT)
  include_directories(${ROOT_INCLUDE_DIR})
  include(${ROOT_USE_FILE})
  add_definitions(-DEUDAQ_LIB_ROOT6)

  # adding command line interface to AidaTLU
  set(EXE_AIDA_TLU_CONTROL AidaTluControl)
  add_executable(${EXE_AIDA_TLU_CONTROL} src/AidaTluControl.cxx)
  target_link_libraries(${EXE_AIDA_TLU_CONTROL} ${EUDAQ_CORE_LIBRARY} ${USER_HARDWARE_LIBRARY} ${ROOT_LIBRARIES})
  list(APPEND INSTALL_TARGETS ${EXE_AIDA_TLU_CONTROL})

  # adding command line interface to AidaTLU
  set(EXE_AIDA_PMTOPTI AidaOptimizeTrigger)
  add_executable(${EXE_AIDA_PMTOPTI} src/AidaOptimizeTrigger.cxx)
  target_link_libraries(${EXE_AIDA_PMTOPTI} ${EUDAQ_CORE_LIBRARY} ${USER_HARDWARE_LIBRARY} ${ROOT_LIBRARIES})
  list(APPEND INSTALL_TARGETS ${EXE_AIDA_PMTOPTI})


  # adding command line interface to AidaTLU
  set(EXE_AIDA_CLI AidaClInterface)
  add_executable(${EXE_AIDA_CLI} src/AidaClInterface.cxx)
  target_link_libraries(${EXE_AIDA_CLI} ${EUDAQ_CORE_LIBRARY} ${USER_HARDWARE_LIBRARY})
  list(APPEND INSTALL_TARGETS ${EXE_AIDA_CLI})
  

else()
  message(STATUS "AIDA TLU stand-alone executables are NOT to be built (USER_BUILD_AIDA_TLU=OFF)")
endif()




install(TARGETS ${INSTALL_TARGETS}
  DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
