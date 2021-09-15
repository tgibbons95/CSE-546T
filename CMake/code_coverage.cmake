function(ENABLE_CODE_COVERAGE TARGET_NAME)
   if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      # using Clang
   elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      # using GCC     
      target_compile_options(${TARGET_NAME} PUBLIC $<$<CONFIG:DEBUG>:-fprofile-arcs> $<$<CONFIG:DEBUG>:-ftest-coverage>)
      target_link_options(${TARGET_NAME} PUBLIC $<$<CONFIG:DEBUG>:-lgcov> $<$<CONFIG:DEBUG>:--coverage>)
   elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      # using Intel C++
   elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
      # using Visual Studio C++
   endif()
endfunction()

function(INIT_CODE_COVERAGE)
   # Include directories
   if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      # using Clang
   elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      # using GCC
      include_directories(/usr/include)
      link_directories(/usr/lib)
   elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      # using Intel C++
   elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
      # using Visual Studio C++
   endif()
endfunction()