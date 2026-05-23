cmake_minimum_required(VERSION 3.10)
project(ekf2 VERSION 0)

# Set this source directory as the tree root for project.
set(ae2f_ProjRoot ${CMAKE_CURRENT_SOURCE_DIR} CACHE STRING "Tree root for project")
set(ae2f_BinRoot ${CMAKE_CURRENT_BINARY_DIR} CACHE STRING "Tree root for project")
set(ae2f_submod  .submod CACHE STRING "Relative path to submodule (subdirectory)")
set(ae2f_CoreScript_Wh ${ae2f_ProjRoot}/${ae2f_submod}/ae2f/CoreScript)

# Execute git to fetch CoreScript when repository does not exist
if(NOT EXISTS ${ae2f_CoreScript_Wh})
	execute_process(
		COMMAND git clone https://talkischeap.work.gd/ae2f/CoreScript
		${ae2f_CoreScript_Wh}
		)
endif()

# add subdirectory
if(NOT TARGET ae2f::CoreScript)
	add_subdirectory(${ae2f_CoreScript_Wh} ${ae2f_BinRoot}/${ae2f_submod}/ae2f/CoreScript)
endif()

ae2f_CoreLibFetch(ae2f Core main)

add_library(Vogen-ekf2 INTERFACE)
add_library(Vogen::ekf2 ALIAS Vogen-ekf2)
target_include_directories(
	Vogen-ekf2 
	INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc/>
	)
target_link_libraries(Vogen-ekf2
	INTERFACE ae2f::Core
	INTERFACE m
	)

ae2f_CoreTestTent(Vogen-ekf2 test)
