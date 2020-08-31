set(sources
  nanoflann.f90
  nanoflann.cpp
  )

add_library(nanoflann_objlib OBJECT ${sources})

include_directories(origin/include)
set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(nanoflann_objlib PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})
set_target_properties(nanoflann_objlib PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

target_include_directories(nanoflann_objlib PUBLIC ${includedir})
