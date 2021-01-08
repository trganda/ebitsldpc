# Attention, this script only support on unix-like system

# define the MATIO_STATIC macro if static build was chosen
if(MATIO_STATIC_LIBRARIES)
    add_definitions(-DMATIO_STATIC)
endif()

# define the list of search paths for headers and libraries
set(FIND_MATIO_PATHS
        ${MATIO_ROOT}
        $ENV{MATIO_ROOT}
        /usr/local/matio
        /usr/local
        /usr
        /opt/hdf5
        /opt/local
        /opt)

# find the matio.h
set(MATIO_HEADER_FOUND FALSE)
find_path(MATIO_INCLUDE_DIR matio.h
        PATH_SUFFIXES include
        PATHS ${FIND_MATIO_PATHS})

# check
if(MATIO_INCLUDE_DIR)
    set(MATIO_HEADER_FOUND TRUE)
endif()

# find the libraries
set(MATIO_FOUND FALSE)
if(MATIO_STATIC_LIBRARIES)
    find_library(MATIO_LIBRARY_STATIC libmatio.a
            PATH_SUFFIXES lib lib64
            PATHS ${FIND_MATIO_PATHS})
    find_library(HDF5_LIBRARY_STATIC libhdf5.a
            PATH_SUFFIXES lib lib64
            PATHS ${FIND_MATIO_PATHS})
    set(MATIO_LIBRARY ${MATIO_LIBRARY_STATIC} ${HDF5_LIBRARY_STATIC})
else()
    find_library(MATIO_LIBRARY_DYNAMIC matio
            PATH_SUFFIXES lib lib64
            PATHS ${FIND_MATIO_PATHS})
    set(MATIO_LIBRARY ${MATIO_LIBRARY_DYNAMIC})
endif()

# check
if(MATIO_LIBRARY AND MATIO_HEADER_FOUND)
    set(MATIO_FOUND TRUE)
endif()

if(MATIO_FOUND)
    message(STATUS "Found matio in ${MATIO_INCLUDE_DIR}")
else()
    message(STATUS "Could NOT found matio")
endif()
