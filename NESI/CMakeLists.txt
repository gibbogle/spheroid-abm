# CMake project file for spheroid-abm

cmake_minimum_required (VERSION 2.6)
project (spheroid-abm)

set (CMAKE_Fortran_COMPILER "ifort")
#set (CMAKE_C_COMPILER "icc")
#set (CMAKE_CXX_COMPILER "icpc")

enable_language (Fortran)
#enable_language (C)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

#set (Fortran_COMPILER_NAME "ifort")

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the NORMA directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message ("Fortran compiler: " ${Fortran_COMPILER_NAME})

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fpp -fopenmp")
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME STREQUAL "gfortran")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

set (SOURCES 
  "/home/mbog002/spheroid-abm/src/bdry_ll.f90" 
  "/home/mbog002/spheroid-abm/src/boundary.f90" 
  "/home/mbog002/spheroid-abm/src/cellstate.f90" 
  "/home/mbog002/spheroid-abm/src/chemokine.f90" 
  "/home/mbog002/spheroid-abm/src/deform.f90"
  "/home/mbog002/spheroid-abm/src/drop.f90"
  "/home/mbog002/spheroid-abm/src/global.f90"
  "/home/mbog002/spheroid-abm/src/nleq1.f90"
  "/home/mbog002/spheroid-abm/src/nleq1e.f90"
  "/home/mbog002/spheroid-abm/src/ode_diffuse.f90"
  "/home/mbog002/spheroid-abm/src/par_zig_mod.f90"
  "/home/mbog002/spheroid-abm/src/rkc.f90" 
  "/home/mbog002/spheroid-abm/src/rkf45s.f90"
  "/home/mbog002/spheroid-abm/src/rksuite.f90"
  "/home/mbog002/spheroid-abm/src/spheroid.f90" 
  "/home/mbog002/spheroid-abm/src/spheroid_main.f90" 
  "/home/mbog002/spheroid-abm/src/Winsock_c.f90" 
  "./clientc.o"
  "./spheroid-abm/tcpc.o"
)

# build executables
#set (NMPROGRAMS "diagstd" "diagrtb" "proj_modes_bin" "pdbmat")
#set (EXECUTABLES "NORMA.exe" ${NMPROGRAMS})
#set (SCRIPTS "gen_pert.sh" "pert_multi_mode.sh")

add_library (spheroid-abm SHARED ${SOURCES})
add_executable (spheroid_main ../src/spheroid_main.f90)
target_link_libraries(spheroid_main spheroid-abm)
set_target_properties(spheroid_main PROPERTIES LINKER_LANGUAGE Fortran)

#foreach (p ${NMPROGRAMS})
#  add_executable (${p} "${p}.f")
#endforeach (p)

# install executables and scripts
#install (TARGETS ${EXECUTABLES} 
#         RUNTIME DESTINATION "bin")
#install (PROGRAMS ${SCRIPTS}
#         DESTINATION "bin") 
