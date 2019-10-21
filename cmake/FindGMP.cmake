# Try to find GMP

# This macro will define the following variables

# GMP_FOUND        - GMP is installed on the system
# GMP_INCLUDE_DIRS - The required include directories for GMP
# GMP_LIBRARIES    - The path to the GMP libraries

include(LibFindMacros)

libfind_pkg_check_modules(GMP_PKGCONF gmp)

find_path(GMP_INCLUDE_DIR
	NAMES gmp.h
	PATHS ${GMP_PKGCONF_INCLUDE_DIRS}
)

find_library(GMP_LIBRARY
	NAMES gmp
	PATHS ${GMP_PKGCONF_LIBRARY_DIRS}
)

set(GMP_PROCESS_INCLUDES GMP_INCLUDE_DIR)
set(GMP_PROCESS_LIBS GMP_LIBRARY)
libfind_process(GMP)
