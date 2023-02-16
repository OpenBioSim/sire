# FindRDKit.cmake
# Placed in the public domain by NextMove Software in 2013
# Try to find RDKit headers and libraries
# Defines:
#
#  RDKIT_FOUND - system has RDKit
#  RDKIT_INCLUDE_DIR - the RDKit include directory
#  RDKIT_LIBRARIES - Link these to use RDKit
#
#  This has been modified to only look in the conda directory
#  for the include and library dirs
#
# References:
#
#  http://nextmovesoftware.com/blog/2013/02/04/looking-for-a-c-cheminformatics-toolkit/
#  https://github.com/timvdm/MolDB/blob/master/cmake/modules/FindRDKit.cmake

include(FindPackageHandleStandardArgs)

if(RDKIT_INCLUDE_DIR AND RDKIT_LIBRARIES)
  # in cache already or user-specified
  find_package_handle_standard_args(RDKit  DEFAULT_MSG
                                    RDKIT_INCLUDE_DIR RDKIT_LIBRARIES)
else()

  if(NOT RDKIT_INCLUDE_DIR)
    find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
              PATHS
              ${CONDA_INCLUDE_DIR}/rdkit
             )
  endif()

  if(NOT RDKIT_LIBRARIES)
    find_library(FILEPARSERS_LIB NAMES FileParsers RDKitFileParsers
                 PATHS
                 ${CONDA_LIBRARY_DIR}

                 #ignore default path, so search starts with above paths
                 NO_DEFAULT_PATH
                 )

    #run with default paths this time
    find_library(FILEPARSERS_LIB NAMES FileParsers RDKitFileParsers)

    if(FILEPARSERS_LIB)
       GET_FILENAME_COMPONENT(RDKIT_LIBRARY_DIR ${FILEPARSERS_LIB} PATH)
       message(STATUS "Found RDKit libraries at ${RDKIT_LIBRARY_DIR}")

      # Note that the order of the following libraries is significant!!
      find_library(SMILESPARSE_LIB NAMES SmilesParse RDKitSmilesParse
                                   HINTS ${RDKIT_LIBRARY_DIR})
      find_library(DEPICTOR_LIB NAMES Depictor RDKitDepictor
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(GRAPHMOL_LIB NAMES GraphMol RDKitGraphMol
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(RDGEOMETRYLIB_LIB NAMES RDGeometryLib RDKitRDGeometryLib
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(RDGENERAL_LIB NAMES RDGeneral RDKitRDGeneral
                                 HINTS ${RDKIT_LIBRARY_DIR})

      find_library(SUBSTRUCTMATCH_LIB NAMES SubstructMatch RDKitSubstructMatch
                                 HINTS ${RDKIT_LIBRARY_DIR})
      find_library(SUBGRAPHS_LIB NAMES Subgraphs RDKitSubgraphs
                                 HINTS ${RDKIT_LIBRARY_DIR})
      find_library(DATASTRUCTS_LIB NAMES DataStructs RDKitDataStructs
                                 HINTS ${RDKIT_LIBRARY_DIR})

      find_library(FORCEFIELD_LIB NAMES ForceField RDKitForceField
                                  HINTS ${RDKIT_LIBRARY_DIR})
      find_library(DISTGEOMHELPERS_LIB NAMES DistGeomHelpers RDKitDistGeomHelpers
                                  HINTS ${RDKIT_LIBRARY_DIR})

      set (RDKIT_LIBRARIES ${FILEPARSERS_LIB} ${SMILESPARSE_LIB}
              ${SUBSTRUCTMATCH_LIB} ${GRAPHMOL_LIB} ${RDGEOMETRYLIB_LIB} ${RDGENERAL_LIB}
              ${SUBGRAPHS_LIB} ${DATASTRUCTS_LIB} ${DEPICTOR_LIB} ${FORCEFIELD_LIB}
              ${DISTGEOMHELPERS_LIB}
              )
    endif()
    if(RDKIT_LIBRARIES)
            message(STATUS "Found RDKit library files at ${RDKIT_LIBRARIES}")
    endif()
  endif()

  find_package_handle_standard_args(RDKit  DEFAULT_MSG
                                  RDKIT_INCLUDE_DIR RDKIT_LIBRARIES)
  mark_as_advanced(RDKIT_INCLUDE_DIR RDKIT_LIBRARIES RDKIT_LIBRARY_DIR)
endif()
