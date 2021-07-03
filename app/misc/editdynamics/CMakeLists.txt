set(projectdir ${PROJECT_SOURCE_DIR})

set(sources-fpp inputdata.F90 parser.F90 editrestart.F90)

set(fypp-opts -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(editrestart ${sources-f90-preproc})

target_link_libraries(editrestart PRIVATE dftbplus)

install(TARGETS editrestart DESTINATION "${CMAKE_INSTALL_BINDIR}")

file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../" projectdir)
