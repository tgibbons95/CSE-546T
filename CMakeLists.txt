cmake_minimum_required(VERSION 3.16)

project(CSE-546T)

enable_testing()

# Include diretories
include_directories(/usr/include)
link_directories(/usr/lib)

add_subdirectory(ExternalSource)

add_subdirectory(BuildingBlocks)
add_subdirectory(Algorithm)
add_subdirectory(Demo)