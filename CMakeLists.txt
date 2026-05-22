cmake_minimum_required(VERSION 3.13)

project(jabber LANGUAGES CXX VERSION 1.0.0)

# Main Jabber compiler options:
option(JABBER_BUILD_PARTICIPANT 
   "Build Jabber preCICE participant app (requires preCICE)." OFF)

option(JABBER_BUILD_VIZ 
   "Build Jabber flowfield visualizer app (requires MFEM)." OFF)

option(JABBER_BUILD_PROFILER 
   "Build Jabber compute profiler app." ON)

option(JABBER_BUILD_OUT 
   "Build Jabber output app." ON)

option(JABBER_BUILD_PLOT 
   "Build Jabber plot app." ON)

option(JABBER_BUILD_PSD
   "Build Jabber PSD app." ON)

option(JABBER_BUILD_DOCS 
   "Build Jabber documentation w/ Doxygen." OFF)

option(JABBER_BUILD_APP_LIB 
   "Build the jabber_app library. Required to build apps." ON)

option(JABBER_BUILD_TESTS 
   "Build unit tests" ON)

option(JABBER_ENABLE_OPENMP
   "Build jabber with OpenMP support" OFF)

option(JABBER_ENABLE_MPI
   "Build jabber_app with MPI support" OFF)

if (JABBER_ENABLE_OPENMP)
   find_package(OpenMP REQUIRED COMPONENTS CXX)
   add_compile_definitions(JABBER_WITH_OPENMP)
endif()

if (JABBER_ENABLE_MPI)
   find_package(MPI REQUIRED COMPONENTS CXX)
   add_compile_definitions(JABBER_WITH_MPI)
endif()

# If building app library
if (JABBER_BUILD_APP_LIB)
   add_compile_definitions(JABBER_WITH_APP)

   # Add TOML11 -- required for app library
   add_subdirectory(externals/toml11)
endif()

# Create the jabber library
add_library(jabber)
set_property(GLOBAL APPEND PROPERTY JABBER_TARGETS_PROPERTY jabber)
target_compile_features(jabber PUBLIC cxx_std_20)
add_subdirectory(jabber)
target_sources(jabber
   PRIVATE 
      ${JABBER_SOURCES}
   PUBLIC
      FILE_SET public_headers TYPE HEADERS
      FILES
         ${JABBER_HEADERS})

# Install jabber
include(GNUInstallDirs)

message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")

# From: https://stackoverflow.com/questions/17511496/how-to-create-a-shared-library-with-cmake
install(TARGETS jabber
   EXPORT jabberTargets
   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   FILE_SET public_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT jabberTargets
   FILE "jabber-targets.cmake"
   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jabber)
   
include(CMakePackageConfigHelpers)
configure_package_config_file("cmake/jabber-config.cmake.in"
   "jabber-config.cmake"
   INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jabber)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/jabber-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES 
   ${CMAKE_CURRENT_BINARY_DIR}/jabber-config.cmake
   ${CMAKE_CURRENT_BINARY_DIR}/jabber-config-version.cmake
   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jabber)

# If building apps
if (JABBER_BUILD_PARTICIPANT 
      OR JABBER_BUILD_VIZ 
      OR JABBER_BUILD_PROFILER
      OR JABBER_BUILD_OUT 
      OR JABBER_BUILD_PLOT
      OR JABBER_BUILD_PSD)

   # All apps use cxxopts
   add_subdirectory(externals/cxxopts)

   # Add app subdirectory
   add_subdirectory(apps)
   
endif()

if (JABBER_BUILD_DOCS)
   add_subdirectory(docs)
endif()

if (JABBER_BUILD_TESTS)
   add_subdirectory(externals/Catch2)
   add_subdirectory(tests)
endif()
