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
   get_cxx_files(MEXMEM_MAIN_HEADER_FILE_PATHS ${SUBMODULES_DIR}/MexMemoryInterfacing/Headers)
   get_cxx_files(MEXMEM_FVT_HEADER_FILE_PATHS ${SUBMODULES_DIR}/MexMemoryInterfacing/Headers/FlatVectTree)
   get_cxx_files(RANDNUMGEN_HEADER_FILE_PATHS ${SUBMODULES_DIR}/RandomNumGen/Headers)
   get_cxx_files(EXTINCURR_HEADER_FILE_PATHS ${SUBMODULES_DIR}/ExternalInputCurrent)
   get_cxx_files(GRID2D_HEADER_FILE_PATHS ${SUBMODULES_DIR}/Grid2D/Headers ${SUBMODULES_DIR}/Grid2D/MexHeaders)
   get_cxx_files(TDNS_HEADER_FILE_PATHS TimeDelNetSim/Headers)
   get_cxx_files(RESPSPIKES_HEADER_FILE_PATHS getResponsibleSpikes/Headers)

   set(COMMON_SOURCE_FILES
         ${MEXMEM_MAIN_HEADER_FILE_PATHS}
         ${MEXMEM_FVT_HEADER_FILE_PATHS}
         ${RANDNUMGEN_HEADER_FILE_PATHS})

   set(TDNS_COMMON_SOURCE_FILES
         ${COMMON_SOURCE_FILES}
         ${EXTINCURR_HEADER_FILE_PATHS}
         ${TDNS_HEADER_FILE_PATHS})

   set(RESPSPIKES_COMMON_SOURCE_FILES
         ${COMMON_SOURCE_FILES}
         ${GRID2D_HEADER_FILE_PATHS}
         ${RESPSPIKES_HEADER_FILE_PATHS})

   set(TDNS_SOURCE_FILES_LIB
         TimeDelNetSim/Source/MexFunctionInterface.cpp 
         TimeDelNetSim/Source/ParallelStructuresII_mex.cpp 
         ${TDNS_COMMON_SOURCE_FILES})
   set(TDNS_SOURCE_FILES_EXE
         TimeDelNetSim/Source/MexFunctionInterfaceExe.cpp
         TimeDelNetSim/Source/ParallelStructuresII_mex.cpp 
         ${TDNS_COMMON_SOURCE_FILES})

   set(RESPSPIKES_SOURCE_FILES_EXE
         getResponsibleSpikes/Source/getRespSpikes_MexInterface.cpp
         getResponsibleSpikes/Source/getRespSpikes_ExeInterface.cpp
         ${RESPSPIKES_COMMON_SOURCE_FILES})
   set(RESPSPIKES_SOURCE_FILES_LIB
         getResponsibleSpikes/Source/getRespSpikes_MexInterface.cpp
         ${RESPSPIKES_COMMON_SOURCE_FILES})

   set(SOURCE_FILES 
         ${TDNS_SOURCE_FILES_LIB}
         ${TDNS_SOURCE_FILES_EXE}
         ${RESPSPIKES_SOURCE_FILES_LIB}
         ${RESPSPIKES_SOURCE_FILES_EXE}
         ${EXAMPLE_DYN_SYS_HEADERS})

   # source group information to generate filters for VS
   source_group(
      "Header Files\\MexMemoryInterfacing" 
      FILES ${MEXMEM_HEADER_FILE_PATHS})
   source_group(
      "Header Files\\RandomNumGen" 
      FILES ${RANDNUMGEN_HEADER_FILE_PATHS})
   source_group(
      "Header Files\\Grid2D"
      FILES ${GRID2D_HEADER_FILE_PATHS})


# assign default build type in case of single config generator 
   if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
      message(STATUS "CMAKE_BUILD_TYPE not defined")
      set(CMAKE_BUILD_TYPE "Release")
   endif()
   message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

# Add Mex_LIB and MEX_EXE Targets for TimeDelNetSim
   custom_add_mex(TimeDelNetSim_mexlib ${TDNS_SOURCE_FILES_LIB})
   enable_tbb(TimeDelNetSim_mexlib)
   custom_add_mex_exec(TimeDelNetSim_mexexe ${TDNS_SOURCE_FILES_EXE})
   enable_tbb(TimeDelNetSim_mexexe)

# Add Mex_LIB and MEX_EXE Targets for getResponsibleSpikes
   custom_add_mex(getResponsibleSpikes_mexlib ${RESPSPIKES_SOURCE_FILES_LIB}
      DOCUMENTATION ${CMAKE_CURRENT_SOURCE_DIR}/getResponsibleSpikes/Docs/MexFunctionDoc.txt)
   enable_tbb(getResponsibleSpikes_mexlib)
   custom_add_mex_exec(getResponsibleSpikes_mexexe ${RESPSPIKES_SOURCE_FILES_EXE})
   enable_tbb(getResponsibleSpikes_mexexe)

# Add Mex_LIB for binarySearch
   custom_add_mex(binarySearch 
      ${CMAKE_SOURCE_DIR}/TimeDelNetSim/MatlabSource/binarySearch/binarySearch.cpp
      DOCUMENTATION ${CMAKE_SOURCE_DIR}/TimeDelNetSim/MatlabSource/binarySearch/binarySearch.txt)

# Add Module Directory File
   configure_file(
      ModulePath.txt.in
      ModulePath.txt)
   
# install instructions
   
   install(DIRECTORY DESTINATION "${CMAKE_BINARY_DIR}/install" CONFIGURATIONS "Release")
   set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")
   install_mex(TimeDelNetSim_mexlib .)
   install_mex_exec(TimeDelNetSim_mexexe .)

   install(FILES ${CMAKE_BINARY_DIR}/ModulePath.txt DESTINATION .)

   install_mex(getResponsibleSpikes_mexlib .)
   install_mex_exec(getResponsibleSpikes_mexexe .)

   install_mex(binarySearch ${CMAKE_SOURCE_DIR}/TimeDelNetSim/MatlabSource)
   
# unset some cached variables so that they are restord to their
# default values the next time cmake is run
   unset(SUBMODULES_DIR CACHE)
