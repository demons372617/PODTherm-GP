# Set PETSc installation path
set(PETSC_DIR "/home/junli/petsc")

# Set the required environment variables
set(ENV{PATH} "${PETSC_DIR}/bin:$ENV{PATH}")
set(ENV{LD_LIBRARY_PATH} "${PETSC_DIR}/lib:$ENV{LD_LIBRARY_PATH}")
set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
set(ENV{C_INCLUDE_PATH} "${PETSC_DIR}/include:$ENV{C_INCLUDE_PATH}")

# Set additional compiler flags if needed
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive")

# Include MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Include PETSc
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDE_DIRS})