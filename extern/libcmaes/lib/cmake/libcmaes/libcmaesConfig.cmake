include (CMakeFindDependencyMacro)
list (APPEND CMAKE_MODULE_PATH /Users/nullptr/Desktop/libcmaes/cmake ${CMAKE_CURRENT_LIST_DIR})

find_dependency (Eigen3 REQUIRED)
if (ON AND TRUE)
  find_dependency (OpenMP REQUIRED)
endif ()

include ("${CMAKE_CURRENT_LIST_DIR}/libcmaesTargets.cmake")
