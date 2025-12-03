#======================================================================================================================
#! CMAKELIST.TXT CORREGIDO - Orden correcto de declaraciones
#======================================================================================================================

# Requerir la versión mínima de CMake.
cmake_minimum_required(VERSION 3.16)

# Definir el nombre del proyecto.
project(barrier VERSION 1.1.0 LANGUAGES CXX)

# Establecer el estándar de C++.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP encontrado - Paralelización habilitada")
endif()

# Establecer el tipo de build por defecto
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "ReleaseFast")
endif()

#======================================================================================================================
#! 1️⃣ CREAR EJECUTABLE PRIMERO (antes de configurarlo)
#======================================================================================================================
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

# Propiedades básicas del target
target_include_directories(barrier PRIVATE include)

#======================================================================================================================
#! 2️⃣ CONFIGURAR FLAGS DE COMPILACIÓN (después de crear el target)
#======================================================================================================================

# Opciones de compilación modernas
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(barrier PRIVATE
        -g -O0 -Wall -Wextra -Wpedantic
        -ggdb3 
        -fno-omit-frame-pointer
        -fno-optimize-sibling-calls
        -fsanitize=address           # Detecta buffer overflows
        -fsanitize=undefined         # Detecta undefined behavior
    )
    target_link_options(barrier PRIVATE 
        -fsanitize=address 
        -fsanitize=undefined
    )
    target_compile_definitions(barrier PRIVATE DEBUG_MODE)

elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    # ✅ OPTIMIZADO: Máximo rendimiento con paralelización
    target_compile_options(barrier PRIVATE
        -O3 -DNDEBUG -Wall -Wextra 
        -march=native -mtune=native  # Optimizar para CPU actual
        -ffast-math                  # Matemática rápida
        -funroll-loops -fomit-frame-pointer
        -flto                        # Link Time Optimization
        -fno-signed-zeros            # Optimizaciones flotantes
        -fno-trapping-math
        -fassociative-math
    )
    target_link_options(barrier PRIVATE
        -flto=auto                   # LTO paralelo
        -Wl,--gc-sections            # Eliminar secciones no usadas
        -Wl,--strip-all              # Strip símbolos debug
    )
    target_compile_definitions(barrier PRIVATE 
        NDEBUG                    # Desactiva assert()
        EIGEN_NO_DEBUG           # Si usas Eigen
        EIGEN_NO_STATIC_ASSERT   # Optimizaciones Eigen
    )

elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    # ✅ COMPROMISO: Optimizado pero con debug info para profiling
    target_compile_options(barrier PRIVATE
        -O3 -g -DNDEBUG -Wall -Wextra
        -funroll-loops
        -flto
    )
    target_link_options(barrier PRIVATE -flto)
    target_compile_definitions(barrier PRIVATE 
        NDEBUG
        EIGEN_NO_DEBUG
        EIGEN_NO_STATIC_ASSERT
    )

elseif(CMAKE_BUILD_TYPE STREQUAL "ReleaseFast")
    # ✅ MÁXIMO RENDIMIENTO: Solo cuando sepas que funciona
    target_compile_options(barrier PRIVATE
        -O3 -DNDEBUG -march=native -ffast-math 
        -funroll-loops -fomit-frame-pointer
        -flto
    )
    target_link_options(barrier PRIVATE -flto=auto -Wl,--gc-sections -Wl,--strip-all)
    target_compile_definitions(barrier PRIVATE NDEBUG)

elseif(CMAKE_BUILD_TYPE STREQUAL "PGOGenerate")
    # 🔍 PASO 1: Generar datos de perfil
    target_compile_options(barrier PRIVATE
        -O3 -DNDEBUG -march=native -ffast-math 
        -funroll-loops -fomit-frame-pointer
        -fprofile-generate           # Instrumentación PGO
    )
    target_link_options(barrier PRIVATE -fprofile-generate)
    target_compile_definitions(barrier PRIVATE NDEBUG)
    message(STATUS "🔍 PGO: Modo generación de perfil activado")

elseif(CMAKE_BUILD_TYPE STREQUAL "PGOUse")
    # ⚡ PASO 2: Usar datos de perfil para optimización máxima
    target_compile_options(barrier PRIVATE
        -O3 -DNDEBUG -march=native -ffast-math 
        -funroll-loops -fomit-frame-pointer
        -flto
        -fprofile-use                # Optimización basada en perfil
        -fprofile-correction         # Corregir inconsistencias menores
    )
    target_link_options(barrier PRIVATE 
        -fprofile-use 
        -flto=auto 
        -Wl,--gc-sections 
        -Wl,--strip-all
    )
    target_compile_definitions(barrier PRIVATE NDEBUG)
    message(STATUS "⚡ PGO: Modo optimización con perfil activado")

else()
    target_compile_options(barrier PRIVATE -O2 -g -Wall -Wextra)
endif()

#======================================================================================================================
#! 3️⃣ BUSCAR Y CONFIGURAR DEPENDENCIAS
#======================================================================================================================

# Buscar dependencias usando métodos modernos
find_package(PkgConfig REQUIRED)

# NetCDF (librería científica)
pkg_check_modules(NETCDF REQUIRED netcdf)
if(NETCDF_FOUND)
    target_include_directories(barrier PRIVATE ${NETCDF_INCLUDE_DIRS})
    target_link_libraries(barrier PRIVATE ${NETCDF_LIBRARIES})
    target_compile_definitions(barrier PRIVATE ${NETCDF_CFLAGS_OTHER})
endif()

# Librerías matemáticas
target_link_libraries(barrier PUBLIC m)

#======================================================================================================================
#! 4️⃣ CONFIGURACIÓN ESPECÍFICA DE COMPILADORES
#======================================================================================================================

# Configuración específica para compiladores
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(barrier PUBLIC -fconcepts-diagnostics-depth=3)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(barrier PUBLIC -fcolor-diagnostics)
endif()

# ⚡ Habilitar OpenMP para paralelización
if(OpenMP_CXX_FOUND)
    target_link_libraries(barrier PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP habilitado y enlazado")
else()
    message(WARNING "OpenMP no encontrado - ejecución secuencial")
endif()

#======================================================================================================================
#! 5️⃣ INFORMACIÓN DE BUILD
#======================================================================================================================

# Información de build
message(STATUS "═══════════════════════════════════════")
message(STATUS "Project: ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "═══════════════════════════════════════")
