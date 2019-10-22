####################################################################################################
# Build Tests
####################################################################################################
project(TransitivityClusteringILPTests)

include_directories("${GTEST_SRC_DIR}/include")

include(cmake/TestExecuter.cmake)
create_test_config_file()

add_library(test_driver main_test.cpp)
COMPILE_FLAGS(test_driver)
target_include_directories(test_driver PUBLIC
	"${GTEST_SRC_DIR}/include"
	${Boost_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
	"${CMAKE_SOURCE_DIR}"
	"${CMAKE_BINARY_DIR}"
)

# Füge einen Test hinzu. Sie können die folgende
# Zeile für neue Tests einfach kopieren und anpassen.
add_gtest(Statistic_test LIBRARIES clustLib)
add_gtest(Distance_test LIBRARIES clustLib)
