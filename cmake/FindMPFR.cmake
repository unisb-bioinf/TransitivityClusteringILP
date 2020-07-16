# Try to find MPFR

# This macro will define the following variables

# MPFR_FOUND        - MPFR is installed on the system
# MPFR_INCLUDE_DIRS - The required include directories for MPFR
# MPFR_LIBRARIES    - The path to the MPFR libraries

include(LibFindMacros)

libfind_pkg_check_modules(MPFR_PKGCONF mpfr)

find_path(MPFR_INCLUDE_DIR
	NAMES mpfr.h
	PATHS ${MPFR_PKGCONF_INCLUDE_DIRS}
)

find_library(MPFR_LIBRARY
	NAMES mpfr
	PATHS ${MPFR_PKGCONF_LIBRARY_DIRS}
)

set(MPFR_PROCESS_INCLUDES MPFR_INCLUDE_DIR)
set(MPFR_PROCESS_LIBS MPFR_LIBRARY)
libfind_process(MPFR)
