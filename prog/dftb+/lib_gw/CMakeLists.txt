set(curdir "lib_gw")

set(sources-fpp
  ${curdir}/gwapprox.F90
  ${curdir}/gw_gammamat.F90
  ${curdir}/gw_input.F90
  ${curdir}/gwxandxc.F90
  ${curdir}/ppmodel.F90
  ${curdir}/randomphase.F90
  ${curdir}/transdens.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
