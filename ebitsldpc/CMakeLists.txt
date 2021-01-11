# CMakeList.txt : CMake project for kmldpc, include source and define
set(THREADS_PREFER_PTHREAD_FLAG ON)
# Will add __FILENAME__ macros for all source files, which is the filename without full find_path
# Courtesy of SO
function(define_file_basename_for_sources targetname)
    get_target_property(source_files "${targetname}" SOURCES)
    foreach(sourcefile ${source_files})
        # Get source file's current list of compile definitions.
        get_property(defs SOURCE "${sourcefile}"
                PROPERTY COMPILE_DEFINITIONS)
        # Add the FILE_BASENAME=filename compile definition to the list.
        get_filename_component(basename "${sourcefile}" NAME)
        list(APPEND defs "__FILENAME__=\"${basename}\"")
        # Set the updated compile definitions on the source file.
        set_property(
                SOURCE "${sourcefile}"
                PROPERTY COMPILE_DEFINITIONS ${defs})
    endforeach()
endfunction()

# include src
include_directories(src, include)
aux_source_directory(src SRC)

# Copy config folder
file(COPY ${CMAKE_SOURCE_DIR}/config/ DESTINATION ${PROJECT_BINARY_DIR}/ebitsldpc FILES_MATCHING PATTERN "*")
# Crate records folder
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/ebitsldpc/records)
# Crate logger folder
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/ebitsldpc/logs)

set(BUILD_STATIC FALSE CACHE STRING "Set this to link external libraries statically")

if(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -g")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
endif()

if(WIN32)
    set(BUILD_STATIC TRUE)
endif(WIN32)

# Set static if BUILD_STATIC is set
if (BUILD_STATIC)
    set(MATIO_STATIC_LIBRARIES TRUE)
    # Link libgcc and libstc++ statically as well
    if(CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++ -lwinpthread")
    endif()
endif()

find_package(Threads REQUIRED)

# include lab code
add_subdirectory("lib")

# Add source to this project's executable.
add_executable (ebitsldpc "ebitsldpc.cc" ${SRC})

target_link_libraries(ebitsldpc PRIVATE lab toml11 Threads::Threads)
# add the binary tree to the search path for include files
target_include_directories(ebitsldpc PUBLIC
        "${PROJECT_SOURCE_DIR}/ebitsldpc/lib/lab"
        "${PROJECT_SOURCE_DIR}/ebitsldpc/lib/toml11"
        )

define_file_basename_for_sources(ebitsldpc)
# TODO: Add tests and install targets if needed.
