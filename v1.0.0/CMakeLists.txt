####################
#  OscProbability  #
#      v1.0.0      #
####################

cmake_minimum_required (VERSION 3.20 FATAL_ERROR)

#Project:
project(OscProbability)
set(maj "1")
set(min "0")
set(ptc "0" )
set(version "${maj}.${min}.${ptc}")

#Colors
string (ASCII 27 esc)
set(N "${esc}[0;00m")
set(R "${esc}[0;91m")
set(B "${esc}[0;94m")
set(G "${esc}[0;92m")
set(Y "${esc}[0;93m")
set(P "${esc}[0;95m")
set(C "${esc}[0;96m")
set(W "${esc}[0;97m")

#Information:
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/src/PREM.info PREM)
file(GLOB srcs ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx)
file(GLOB hdrs ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.h)
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/INFO info)
message("\n${G}:: Standard Information ::${N}")
message("\n${C}${info}${N}")
message("\n${G}:: Current PREM Information ::${N}")
message("\n${C}${PREM}${N}")

#Include Eigen3:
find_package (Eigen3 3.3.7 REQUIRED NO_MODULE)
if (Eigen3_FOUND)
    set (eigen3_inc "${EIGEN_INCLUDE_DIR}")
if (eigen3_inc)
message (STATUS "Found Eigen3: ${eigen3_inc}")
else()
message (FATAL_ERROR "${Y}\n EIGEN3 is ${R}NOT${Y} set, Try with \n-DEIGEN_INCLUDE_DIR=\$(path-to-'${P}eigen3-include-dir${Y}') ${N}\n")
endif(eigen3_inc)
endif(Eigen3_FOUND)

#Include GSL:
find_package(GSL 2.4)
if (GSL_FOUND)
    set (gsl_inc   ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})
    set (gsl_flags ${GSL_CONFIG_EXECUTABLE})
    set (gsl_libs  ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARY})
else ()
message (FATAL_ERROR "${Y}\n GSL is ${R}NOT${Y} set, Try with \n-DGSL_ROOT_DIR=\$(path-to-'${P}install-dir${Y}') ${N}\n")
endif(GSL_FOUND)

#Include OpenMP:
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (omp_flags ${OpenMP_CXX_FLAGS})
    set (omp_inc   ${OpenMP_CXX_INCLUDE_DIRS})
    set (omp_libs  ${OpenMP_CXX_LIBRARIES})
endif(OPENMP_FOUND)

#Include BLAS:
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    set (blas_flags ${BLAS_LINKER_FLAGS})
    set (blas_libs  ${BLAS_LIBRARIES})
endif(BLAS_FOUND)

#Include LAPACK:
find_package(LAPACK REQUIRED)
if (LAPACK_FOUND) 
    set (lpk_flags ${LAPACK_LINKER_FLAGS})
    set (lpk_libs  ${LAPACK_LIBRARIES})
endif(LAPACK_FOUND)

#Include ATLAS:
find_library(LIB_ATLAS NAMES atlas)
set(Checklib LIB_ATLAS)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Atlas DEFAULT_MSG ${Checklib})
if(ATLAS_FOUND)
  set(atlas_libs ${LIB_ATLAS})
endif(ATLAS_FOUND)

#Include ROOT (Optional)
option (ROOTSYS OFF)
if (ROOTSYS)
find_package(ROOT 6.18)
if (ROOT_FOUND)
    set (root_inc   ${ROOT_INCLUDE_DIRS})
    set (root_libs  ${ROOT_LIBRARIES})
    set (root_flags ${ROOT_CXX_FLAGS})
endif (ROOT_FOUND)
else  ()
message ("${Y}To Turn ON CERN ROOT ( >= 6.18); -DROOTSYS=ON${N}\n")
endif (ROOTSYS)
unset (ROOTSYS CACHE)

#Include FLAGS:
add_compile_options(-Wall -Wextra -Wpedantic -Werror)
set(CMAKE_CXX_FLAGS "-std=c++11 -msse4.2 -ffast-math -march=native -O3 -DNDEBUG -DEIGEN_USE_BLAS")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${omp_flags} ${lpk_flags} ${root_flags}")

#Include LIBRARIES:
set (lib "OscProbability${version}")
set (libs ${gsl_libs} ${omp_libs} ${lpk_libs} ${blas_libs} ${root_libs} -L${atlas_libs} pthread)
include_directories(${gsl_inc} ${eigen3_inc} ${omp_inc} ${root_inc} ${CMAKE_CURRENT_SOURCE_DIR}/inc)

#STATIC LIBRARY:
add_library(${lib} ${srcs} ${hdrs})

#TARGET CALL:
set (Target ${SET_TARGET})
set (extn ".cxx")
add_executable(${Target} ${Target}${extn} ${hdrs})
target_link_libraries(${Target} PRIVATE ${lib} ${libs})






