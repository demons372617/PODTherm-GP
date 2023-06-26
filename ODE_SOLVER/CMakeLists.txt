cmake_minimum_required(VERSION 3.5)
project(Calculate_CU)

add_compile_options("-g")
set(pugixml_DIR /home/junli/pugixml)
find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)
#find_package(pugixml REQUIRED)
include_directories(${MPI_INCLUDE_PATH} ${pugixml_INCLUDE_DIRS})

# Add the options to download SuperLU and MUMPS
# set(PETSC_USE_SUPERLU ON CACHE BOOL "" FORCE)
# set(PETSC_USE_MUMPS ON CACHE BOOL "" FORCE)

include_directories(${MPI_INCLUDE_PATH})

set(XML_FILE ../POD_Para.xml)
# Set the path to your PETSc installation
set(PETSC_PATH "/home/junli/petsc")

# Set the PKG_CONFIG_PATH to include the PETSc installation
set(ENV{PKG_CONFIG_PATH} "${PETSC_PATH}/lib/pkgconfig")

#pkg_check_modules(PETSC_PKG REQUIRED IMPORTED_TARGET PETSc)
#find_package(PETSC)
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED) 

set(SRC_DIR src/)
set(MAIN_CPP ${SRC_DIR}main.cpp)
set(SOURCES
    ${SRC_DIR}POD_ODE_Solver.hpp
    ${SRC_DIR}POD_ODE_Solver.cpp
    ${SRC_DIR}Matrix.cpp
    ${SRC_DIR}Matrix.hpp
)

set(ALL_SOURCES
    ${MAIN_CPP}
    ${SOURCES}
)

add_compile_options("-lpthread -fpermissive")
set(GCC_CXX_FLAGS ${GCC_CXX_FLAGS} "-std=c++11 -lpthread -fpermissive")

### set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -Wall -Wextra -Werror")
### set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --with-debugging=yes")

### MAIN PROGRAMS #########################################

add_executable(${PROJECT_NAME} ${ALL_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${SRC_DIR}>
    ${MPI_INCLUDE_PATH}
)

list(APPEND CUST_PETSC_INCLUDE_DIRS "/home/junli/petsc/include")
list(APPEND CUST_PETSC_INCLUDE_DIRS "/home/junli/petsc/arch-linux-c-debug/include")
add_library(cust_petsc SHARED IMPORTED) # or STATIC instead of SHARED
set_target_properties(cust_petsc PROPERTIES
  IMPORTED_LOCATION "/home/junli/petsc/arch-linux-c-debug/lib/libpetsc.so.3.19.2"
  INTERFACE_INCLUDE_DIRECTORIES "${CUST_PETSC_INCLUDE_DIRS}"
)

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${SRC_DIR}>
    ${MPI_INCLUDE_PATH}
    ${pugixml_INCLUDE_DIRS}
)

#target_link_libraries(${PROJECT_NAME} PRIVATE 
#    Threads::Threads
#    ${MPI_LIBRARIES}
#   #PkgConfig::PETSC_PKG
#    ${pugixml_LIBRARIES}
#    
#    #PkgConfig::SLEPC_PKG
#)

target_link_libraries(${PROJECT_NAME} PRIVATE 
    Threads::Threads
    ${MPI_LIBRARIES}
    cust_petsc
    pugixml
    #PkgConfig::SLEPC_PKG
)
