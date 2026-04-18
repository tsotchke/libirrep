# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.16)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" LIBIRREP_VERSION_STR)
string(REGEX REPLACE "-.*$" "" LIBIRREP_VERSION_NUM "${LIBIRREP_VERSION_STR}")

project(libirrep
    VERSION     ${LIBIRREP_VERSION_NUM}
    DESCRIPTION "SO(3)/SU(2)/O(3) irreducible-representation machinery"
    LANGUAGES   C
)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(LIBIRREP_BUILD_TESTS      "Build libirrep tests"      ON)
option(LIBIRREP_BUILD_BENCHMARKS "Build libirrep benchmarks" OFF)
option(LIBIRREP_BUILD_EXAMPLES   "Build libirrep examples"   ON)
option(LIBIRREP_INSTALL          "Install libirrep"          ON)

# ---- sources ------------------------------------------------------------
file(GLOB LIBIRREP_SOURCES CONFIGURE_DEPENDS "src/*.c")

add_library(libirrep ${LIBIRREP_SOURCES})
add_library(libirrep::libirrep ALIAS libirrep)

target_include_directories(libirrep
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_definitions(libirrep
    PRIVATE
        IRREP_BUILDING_DLL
)

target_compile_options(libirrep PRIVATE
    -Wall -Wextra -Wpedantic -Wno-unused-parameter
    -fno-math-errno
)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    target_compile_definitions(libirrep PRIVATE USE_NEON_IF_AVAILABLE)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    target_compile_definitions(libirrep PRIVATE
        USE_SSE42_IF_AVAILABLE
        USE_AVX2_IF_AVAILABLE
    )
endif()

set_target_properties(libirrep PROPERTIES
    VERSION       ${PROJECT_VERSION}
    SOVERSION     ${PROJECT_VERSION_MAJOR}
)

if(NOT WIN32)
    target_link_libraries(libirrep PUBLIC m)
endif()

# ---- tests --------------------------------------------------------------
if(LIBIRREP_BUILD_TESTS)
    enable_testing()
    file(GLOB LIBIRREP_TEST_SOURCES CONFIGURE_DEPENDS "tests/test_*.c")
    foreach(src ${LIBIRREP_TEST_SOURCES})
        get_filename_component(name ${src} NAME_WE)
        add_executable(${name} ${src})
        target_link_libraries(${name} PRIVATE libirrep)
        target_include_directories(${name} PRIVATE tests)
        add_test(NAME ${name} COMMAND ${name})
    endforeach()
endif()

# ---- benchmarks ---------------------------------------------------------
if(LIBIRREP_BUILD_BENCHMARKS)
    file(GLOB LIBIRREP_BENCH_SOURCES CONFIGURE_DEPENDS "benchmarks/bench_*.c")
    foreach(src ${LIBIRREP_BENCH_SOURCES})
        get_filename_component(name ${src} NAME_WE)
        add_executable(${name} ${src})
        target_link_libraries(${name} PRIVATE libirrep)
        target_include_directories(${name} PRIVATE benchmarks)
    endforeach()
endif()

# ---- examples -----------------------------------------------------------
if(LIBIRREP_BUILD_EXAMPLES)
    file(GLOB LIBIRREP_EXAMPLE_SOURCES CONFIGURE_DEPENDS "examples/*.c")
    foreach(src ${LIBIRREP_EXAMPLE_SOURCES})
        get_filename_component(name ${src} NAME_WE)
        add_executable(${name} ${src})
        target_link_libraries(${name} PRIVATE libirrep)
    endforeach()
endif()

# ---- install ------------------------------------------------------------
if(LIBIRREP_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    install(TARGETS libirrep
        EXPORT libirrepTargets
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
    install(DIRECTORY include/irrep
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h"
    )

    configure_file(
        scripts/libirrep.pc.in
        ${CMAKE_CURRENT_BINARY_DIR}/libirrep.pc
        @ONLY
    )
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libirrep.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    )
endif()
