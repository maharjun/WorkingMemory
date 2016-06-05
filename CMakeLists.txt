cmake_minimum_required(VERSION 3.5)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")
include(CppGlob)
include(CustomMexFunctions)
include(TBBFunctions)
project(TimeDelNetSim)

FIND_PACKAGE(Matlab REQUIRED)
FIND_PACKAGE(TBB REQUIRED)

message(STATUS "C Compiler   = ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler = ${CMAKE_CXX_COMPILER}")

# define submodules directory and adding into include path
   if(NOT DEFINED ENV{SUBMODULES_DIR} AND NOT DEFINED SUBMODULES_DIR)
      set(SUBMODULES_DIR ${CMAKE_CURRENT_SOURCE_DIR})
   elseif(NOT DEFINED SUBMODULES_DIR)
      set(SUBMODULES_DIR $ENV{SUBMODULES_DIR})
   endif()
   message(STATUS "SUBMODULES_DIR = ${SUBMODULES_DIR}")
   include_directories(BEFORE ${SUBMODULES_DIR})

# add TimeDelNetSim into include path (for ExternalInputCurrent)
   include_directories(BEFORE ${CMAKE_SOURCE_DIR}/TimeDelNetSim)

# include the relevant source and header files
   get_cxx_files(MEXMEM_HEADER_FILE_PATHS ${SUBMODULES_DIR}/MexMemoryInterfacing/Headers)
   get_cxx_files(RANDNUMGEN_HEADER_FILE_PATHS ${SUBMODULES_DIR}/RandomNumGen/Headers)
   get_cxx_files(EXTINCURR_HEADER_FILE_PATHS ${SUBMODULES_DIR}/ExternalInputCurrent)
   get_cxx_files(HEADER_FILE_PATHS Headers)

   set(COMMON_SOURCE_FILES
         ${MEXMEM_HEADER_FILE_PATHS}
         ${RANDNUMGEN_HEADER_FILE_PATHS}
         ${EXTINCURR_HEADER_FILE_PATHS}
         ${HEADER_FILE_PATHS})

   set(SOURCE_FILES_LIB
         TimeDelNetSim/Source/MexFunctionInterface.cpp 
         TimeDelNetSim/Source/ParallelStructuresII_mex.cpp 
         ${COMMON_SOURCE_FILES})
   set(SOURCE_FILES_EXE
         TimeDelNetSim/Source/MexFunctionInterfaceExe.cpp
         TimeDelNetSim/Source/ParallelStructuresII_mex.cpp 
         ${COMMON_SOURCE_FILES})

   set(SOURCE_FILES 
         ${SOURCE_FILES_LIB}
         ${SOURCE_FILES_EXE}
         ${EXAMPLE_DYN_SYS_HEADERS})

   # source group information to generate filters for VS
   source_group(
      "Header Files\\MexMemoryInterfacing" 
      FILES ${MEXMEM_HEADER_FILE_PATHS})
   source_group(
      "Header Files\\RandomNumGen" 
      FILES ${RANDNUMGEN_HEADER_FILE_PATHS})

# assign default build type in case of single config generator 
   if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
      message(STATUS "CMAKE_BUILD_TYPE not defined")
      set(CMAKE_BUILD_TYPE "Release")
   endif()
   message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

# Add Mex_LIB and MEX_EXE Targets
   custom_add_mex(TimeDelNetSim_mexlib ${SOURCE_FILES_LIB})
   enable_tbb(TimeDelNetSim_mexlib)
   custom_add_mex_exec(TimeDelNetSim_mexexe ${SOURCE_FILES_EXE})
   enable_tbb(TimeDelNetSim_mexexe)

# install instructions
   
   install(DIRECTORY DESTINATION "${CMAKE_BINARY_DIR}/install" CONFIGURATIONS "Release")
   set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")
   install_mex(TimeDelNetSim_mexlib .)
   install_mex_exec(TimeDelNetSim_mexexe .)
   
# unset some cached variables so that they are restord to their
# default values the next time cmake is run
   unset(SUBMODULES_DIR CACHE)
