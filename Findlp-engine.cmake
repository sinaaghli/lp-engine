# - Check for the presence of lp-engine
#
# The following variables are set when <PACKAGE> is found:
#  HAVE_lp-engine       = Set to true, if all components of lp-engine
#                          have been found.
#  lp-engine_INCLUDES   = Include path for the header files of lp-engine
#  lp-engine_LIBRARIES  = Link these to use lp-engine

## -----------------------------------------------------------------------------
## Check for the header files

find_path (lp-engine_INCLUDES lpengine.h
  PATHS /usr/local/include /usr/include /sw/include
  PATH_SUFFIXES lp-engine
  )

## -----------------------------------------------------------------------------
## Check for the library

find_library (lp-engine_LIBRARIES lp-engine
  PATHS /usr/local/lib /usr/lib /lib /sw/lib
  )

## -----------------------------------------------------------------------------
## Actions taken when all components have been found

if (lp-engine_INCLUDES AND lp-engine_LIBRARIES)
  set (HAVE_lp-engine TRUE)
else (lp-engine_INCLUDES AND lp-engine_LIBRARIES)
  if (NOT lp-engine_FIND_QUIETLY)
    if (NOT lp-engine_INCLUDES)
      message (STATUS "Unable to find lp-engine header files!")
    endif (NOT lp-engine_INCLUDES)
    if (NOT lp-engine_LIBRARIES)
      message (STATUS "Unable to find lp-engine library files!")
    endif (NOT lp-engine_LIBRARIES)
  endif (NOT lp-engine_FIND_QUIETLY)
endif (lp-engine_INCLUDES AND lp-engine_LIBRARIES)

if (HAVE_lp-engine)
  if (NOT lp-engine_FIND_QUIETLY)
    message (STATUS "Found components for lp-engine")
    message (STATUS "lp-engine_INCLUDES = ${lp-engine_INCLUDES}")
    message (STATUS "lp-engine_LIBRARIES     = ${lp-engine_LIBRARIES}")
  endif (NOT lp-engine_FIND_QUIETLY)
else (HAVE_lp-engine)
  if (lp-engine_FIND_REQUIRED)
    message (FATAL_ERROR "Could not find lp-engine!")
  endif (lp-engine_FIND_REQUIRED)
endif (HAVE_lp-engine)

mark_as_advanced (
  HAVE_lp-engine
  lp-engine_LIBRARIES
  lpengine_INCLUDES
  )
