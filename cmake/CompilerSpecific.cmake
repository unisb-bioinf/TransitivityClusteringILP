####################################################################################################
# Compile flags
####################################################################################################

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.4.0")
		message(FATAL_ERROR "netICS requires a clang version >= 3.4")
	endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0")
		message(FATAL_ERROR "netICS requires a GCC version >= 4.8.0")
	endif()

	# This "fixes" some annoying warnings with Eigen3. We can probably
	# Remove this when newer Eigen versions are available.
	if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "6.0.0")
		LIST(APPEND CXX_FLAGS "-Wno-ignored-attributes")
	endif()
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	LIST(APPEND CXX_FLAGS "-fvisibility=hidden")
	LIST(APPEND CXX_FLAGS "-pedantic")
	LIST(APPEND CXX_FLAGS "-Wall")
	LIST(APPEND CXX_FLAGS "-Wextra")
	LIST(APPEND CXX_FLAGS "-std=c++14")

	if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
		list(APPEND CXX_FLAGS "-fno-omit-frame-pointer")
	endif()

	if(CMAKE_BUILD_TYPE STREQUAL "Release")
		LIST(APPEND CXX_FLAGS "-Ofast")
	endif()

	set(LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
else()
	message(STATUS "At the moment, GeneTrail2 supports only Clang and the GNU Compiler Collection.")
	return()
endif()

## Enable C++11 mode
if    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	LIST(APPEND CXX_FLAGS "-Wno-deprecated-register")
endif()

set(CMAKE_MODULE_LINKER_FLAGS ${LINKER_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS    ${LINKER_FLAGS})
set(CMAKE_SHARED_LINKER_FLAGS ${LINKER_FLAGS})

function(COMPILE_FLAGS target)
	target_compile_options(${target} PUBLIC ${CXX_FLAGS})
endfunction()
