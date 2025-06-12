if (NOT DEFINED PYTHON_EXECUTABLE)
  # we will just use the python that comes with anaconda
  if (MSYS)
      set (PYTHON_EXECUTABLE "${ANACONDA_BASE}/python" )
  elseif (MSVC)
      set (PYTHON_EXECUTABLE "${ANACONDA_BASE}/python.exe" )
  else()
      set (PYTHON_EXECUTABLE "${ANACONDA_BASE}/bin/python3" )
  endif()
endif()

find_package( Python REQUIRED COMPONENTS Interpreter Development )

set( PYTHON_VERSION "${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}" )

message( STATUS "Using python ${PYTHON_VERSION} from ${PYTHON_EXECUTABLE}" )

# Set the require Python variables.
set ( PYTHON_SITE_DIR "${Python_SITELIB}" )
set ( PYTHON_INCLUDE_DIR "${Python_INCLUDE_DIRS}" )
set ( PYTHON_LIBRARIES "${Python_LIBRARIES}" )

if (NOT PYTHON_LIBRARIES)
  message( FATAL_ERROR "Where is the python library that comes with anaconda? "
                        "It cannot be found. Please check that your anaconda "
                        "installation is complete." )
endif()

message( STATUS "Using anaconda/miniconda python in ${PYTHON_LIBRARIES} | ${PYTHON_INCLUDE_DIR}" )
message( STATUS "Python modules will be installed to ${PYTHON_SITE_DIR}" )

set( SIRE_FOUND_PYTHON TRUE )
