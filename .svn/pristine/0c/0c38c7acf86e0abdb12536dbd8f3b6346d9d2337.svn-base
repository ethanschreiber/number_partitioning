cmake_minimum_required(VERSION 2.8.8)
project(Partitioning)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin/")

# Generic build options
set(CMAKE_C_FLAGS                "-Wall -Wextra -pedantic -Wno-long-long")
set(CMAKE_CXX_FLAGS              "-std=c++0x -Wall -Wextra -pedantic -Wno-long-long")
set(CMAKE_EXE_LINKER_FLAGS       "")

# Variant-specific build options
set(CMAKE_C_FLAGS_DEBUG          "-g -O0 -DDEBUG -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_DEBUG        "-g -O0 -DDEBUG -D_GLIBCXX_DEBUG")
set(CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -DNDEBUG")



set(BOOST_PROGRAM_OPTIONS boost_program_options-mt)

add_subdirectory(ss)
add_subdirectory(partition)
add_subdirectory(pack)

############# LIBRARY ################

# Monolithic library of all shared C++ files
#add_library(SHARED_FILES PDBIndexer.cpp EdgePDB.cpp RubiksUtils.cpp FileQueue.cpp)


############# 2-Way Partition ############

	set(partition_source_files
		PartitioningMain.cpp
		MainUtils.cpp
		Utils.cpp
		)
add_executable(Partition ${partition_source_files})
target_link_libraries(Partition lib_ss lib_partition lib_pack ${BOOST_PROGRAM_OPTIONS})

