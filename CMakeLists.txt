# Requerir la versión mínima de CMake.
cmake_minimum_required(VERSION 3.16)

# Definir el nombre del proyecto.
project(barrier VERSION 1.0.0 LANGUAGES CXX)

# Establecer el estándar de C++.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Establecer el tipo de build por defecto
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Opciones de compilación modernas
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -ffast-math")
else()
    # Para cualquier otro build type (como "RelWithDebInfo")
    set(CMAKE_CXX_FLAGS "-O2 -g -Wall -Wextra")
endif()

# Crear ejecutable
add_executable(barrier 
    src/main.cpp 
    src/simulation.cpp 
    src/screen_presenter.cpp 
    src/data_reader.cpp
    src/utils.cpp 
    src/error_handling.cpp 
    src/cross_section.cpp 
    src/hydrograph.cpp 
    src/data_writer.cpp
)

# Propiedades del target
target_include_directories(barrier PRIVATE include)

# Buscar dependencias usando métodos modernos
find_package(PkgConfig REQUIRED)

# NetCDF (librería científica)
pkg_check_modules(NETCDF REQUIRED netcdf)
if(NETCDF_FOUND)
    target_include_directories(barrier PRIVATE ${NETCDF_INCLUDE_DIRS})
    target_link_libraries(barrier PRIVATE ${NETCDF_LIBRARIES})
    target_compile_definitions(barrier PRIVATE ${NETCDF_CFLAGS_OTHER})
endif()

# ✅ ELIMINADO: Todo el código de OpenMP
# Definir que OpenMP está deshabilitado
target_compile_definitions(barrier PRIVATE NO_OPENMP)

# ✅ ELIMINADO: Threads (ya no necesarios sin paralelización)
# find_package(Threads REQUIRED)
# target_link_libraries(barrier PUBLIC Threads::Threads)

# Librerías matemáticas
target_link_libraries(barrier PUBLIC m)

# Opciones adicionales para mejor debugging
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(barrier PRIVATE DEBUG_MODE)
    target_compile_options(barrier PRIVATE 
        -ggdb3 
        -fno-omit-frame-pointer
        -fno-optimize-sibling-calls
    )
endif()

# Configuración específica para compiladores
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(barrier PUBLIC -fconcepts-diagnostics-depth=3)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(barrier PUBLIC -fcolor-diagnostics)
endif()

# Información de build
message(STATUS "═══════════════════════════════════════")
message(STATUS "Project: ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "═══════════════════════════════════════")
