set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/bulkpot.F90
  ${curdir}/fancybc.F90
  ${curdir}/gallocation.F90
  ${curdir}/gewald.F90
  ${curdir}/mpi_poisson.F90
  ${curdir}/parameters.F90
  ${curdir}/parcheck.F90
  ${curdir}/poisson.F90
  ${curdir}/structure.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
