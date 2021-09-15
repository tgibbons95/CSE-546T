project(BuildingBlocks)

add_subdirectory(src)
add_subdirectory(test)

install(DIRECTORY "include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")