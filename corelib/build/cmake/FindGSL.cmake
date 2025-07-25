#
# Try to find gnu scientific library GSL
# (see http://www.gnu.org/software/gsl/)
# Once run this will define:
#
# GSL_FOUND       = system has GSL lib
#
# GSL_LIBRARIES   = full path to the libraries
#    on Unix/Linux with additional linker flags from "gsl-config --libs"
#
# CMAKE_GSL_CXX_FLAGS  = Unix compiler flags for GSL, essentially "`gsl-config --cxxflags`"
#
# GSL_INCLUDE_DIR      = where to find headers
#
# GSL_LINK_DIRECTORIES = link directories, useful for rpath on Unix
# GSL_EXE_LINKER_FLAGS = rpath on Unix
#
# Felix Woelk 07/2004
# minor corrections Jan Woetzel
#
# Developed further by Christopher Woods
#
# www.mip.informatik.uni-kiel.de
# --------------------------------

IF (NOT WIN32)
  SET(GSL_CONFIG_PREFER_PATH
  "$ENV{GSL_DIR}/bin"
  "$ENV{GSL_DIR}"
  "$ENV{GSL_HOME}/bin"
  "$ENV{GSL_HOME}"
  CACHE STRING "preferred path to GSL (gsl-config)")
  FIND_PROGRAM(GSL_CONFIG gsl-config
    ${GSL_CONFIG_PREFER_PATH}
    /usr/bin/
    ENV PATH
    )
  #MESSAGE("DBG GSL_CONFIG ${GSL_CONFIG}")

  IF (GSL_CONFIG)
    # set CXXFLAGS to be fed into CXX_FLAGS by the user:
    EXECUTE_PROCESS( COMMAND ${GSL_CONFIG} --cflags
                  OUTPUT_VARIABLE GSL_CXXFLAGS
                  OUTPUT_STRIP_TRAILING_WHITESPACE )

    # set INCLUDE_DIRS to prefix+include
    EXECUTE_PROCESS( COMMAND ${GSL_CONFIG} --prefix
      OUTPUT_VARIABLE GSL_PREFIX)
    SET(GSL_INCLUDE_DIR ${GSL_PREFIX}/include CACHE STRING INTERNAL)

    # set link libraries and link flags
    EXECUTE_PROCESS( COMMAND ${GSL_CONFIG}
                  --libs-without-cblas
                  OUTPUT_VARIABLE GSL_LIBRARIES 
                  OUTPUT_STRIP_TRAILING_WHITESPACE )

    # extract link dirs for rpath
    EXECUTE_PROCESS( COMMAND ${GSL_CONFIG} --libs
      OUTPUT_VARIABLE GSL_CONFIG_LIBS 
      OUTPUT_STRIP_TRAILING_WHITESPACE )

    # split off the link dirs (for rpath)
    # use regular expression to match wildcard equivalent "-L*<endchar>"
    # with <endchar> is a space or a semicolon
    STRING(REGEX MATCHALL "[-][L]([^ ;])+"
      GSL_LINK_DIRECTORIES_WITH_PREFIX
      "${GSL_CONFIG_LIBS}" )
      #      MESSAGE("DBG  GSL_LINK_DIRECTORIES_WITH_PREFIX=${GSL_LINK_DIRECTORIES_WITH_PREFIX}")

    # remove prefix -L because we need the pure directory for LINK_DIRECTORIES

    IF (GSL_LINK_DIRECTORIES_WITH_PREFIX)
      STRING(REGEX REPLACE "[-][L]" "" GSL_LINK_DIRECTORIES ${GSL_LINK_DIRECTORIES_WITH_PREFIX} )
    ENDIF (GSL_LINK_DIRECTORIES_WITH_PREFIX)
    SET(GSL_EXE_LINKER_FLAGS "-Wl,-rpath,${GSL_LINK_DIRECTORIES}" CACHE STRING INTERNAL)
    #      MESSAGE("DBG  GSL_LINK_DIRECTORIES=${GSL_LINK_DIRECTORIES}")
    #      MESSAGE("DBG  GSL_EXE_LINKER_FLAGS=${GSL_EXE_LINKER_FLAGS}")

    #      ADD_DEFINITIONS("-DHAVE_GSL")
    #      SET(GSL_DEFINITIONS "-DHAVE_GSL")
    MARK_AS_ADVANCED(
      GSL_CXX_FLAGS
      GSL_INCLUDE_DIR
      GSL_LIBRARIES
      GSL_LINK_DIRECTORIES
      GSL_DEFINITIONS
    )
    MESSAGE(STATUS "Using GSL from ${GSL_PREFIX}")

  ELSE(GSL_CONFIG)
    MESSAGE("FindGSL.cmake: gsl-config not found. Please set it manually. GSL_CONFIG=${GSL_CONFIG}")
  ENDIF(GSL_CONFIG)
ENDIF(NOT WIN32)

IF((NOT GSL_LIBRARIES) AND MSVC AND DEFINED SIRE_APP)
  find_library( GSL_LIBRARY "gsl.lib" PATHS "${SIRE_APP}/Library/lib" NO_DEFAULT_PATH )
  find_library( GSLCBLAS_LIBRARY "gslcblas.lib" PATHS "${SIRE_APP}/Library/lib" NO_DEFAULT_PATH )
  set ( GSL_LIBRARIES "${GSL_LIBRARY};${GSLCBLAS_LIBRARY}" )
  get_filename_component(GSL_INCLUDE_DIR "${GSL_LIBRARY}" DIRECTORY)
  set ( GSL_INCLUDE_DIR "${GSL_INCLUDE_DIR}/include" )
ENDIF()

IF(GSL_LIBRARIES)
  IF(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)

    SET(GSL_FOUND 1)

  ENDIF(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)
ENDIF(GSL_LIBRARIES)
