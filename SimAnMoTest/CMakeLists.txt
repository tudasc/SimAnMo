target_sources(SimAnMoLibraryTester
	PRIVATE 
		${CMAKE_CURRENT_LIST_DIR}/main.cpp
)

get_filename_component(PARENT_DIR ../ ABSOLUTE)

target_include_directories(SimAnMo PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           ${PARENT_DIR}/SimulatedAnnealingExtraP
 )


message(STATUS "include dir: ${PARENT_DIR}/SimulatedAnnealingExtraP")
