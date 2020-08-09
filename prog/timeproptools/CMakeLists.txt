set(sources-fpp timepropedit.F90 inputedit.F90)

list(APPEND FYPP_FLAGS -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(timepropedit ${sources-f90-preproc})

if(WITH_MPI)
  target_include_directories(timepropedit PRIVATE ${MPI_Fortran_MODULE_DIR})
  target_link_libraries(timepropedit PRIVATE ${MPI_Fortran_LIBRARIES})
endif()

target_link_libraries(timepropedit PRIVATE dftbplus)

install(TARGETS timepropedit DESTINATION ${INSTALL_BIN_DIR})
