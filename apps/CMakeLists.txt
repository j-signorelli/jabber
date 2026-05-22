set(JABBER_APP_TARGETS "")

# Build participant app
if (JABBER_BUILD_PARTICIPANT)
   # Find preCICE
   find_package(PRECICE 3.1.2 REQUIRED CONFIG NAMES precice)
   message(STATUS "Found preCICE config: TRUE")

   add_executable(jabber_participant jabber_participant.cpp)
   target_link_libraries(jabber_participant PRIVATE 
                                 jabber
                                 cxxopts::cxxopts
                                 precice::precice)
   list(APPEND JABBER_APP_TARGETS jabber_participant)
endif()

# Build Jabber acoustic field visualizer app
if (JABBER_BUILD_VIZ)

   # Find MFEM
   find_package(MFEM 4.8 REQUIRED CONFIG)
   message(STATUS "Found MFEM config: TRUE")

   add_executable(jabber_viz jabber_viz.cpp)
   target_link_libraries(jabber_viz PRIVATE 
                                 jabber
                                 cxxopts::cxxopts
                                 mfem)
   list(APPEND JABBER_APP_TARGETS jabber_viz)
endif()


# Build Jabber profiler app
if (JABBER_BUILD_PROFILER)

   add_executable(jabber_profile jabber_profile.cpp)
   target_link_libraries(jabber_profile PRIVATE 
                                       jabber
                                       cxxopts::cxxopts)   
   list(APPEND JABBER_APP_TARGETS jabber_profile)
endif()


# Build Jabber out app
if (JABBER_BUILD_OUT)

   add_executable(jabber_out jabber_out.cpp)
   target_link_libraries(jabber_out PRIVATE 
                                       jabber
                                       cxxopts::cxxopts)
   list(APPEND JABBER_APP_TARGETS jabber_out)
endif()

# Build Jabber plot app
if (JABBER_BUILD_PLOT)

   add_executable(jabber_plot jabber_plot.cpp)
   target_link_libraries(jabber_plot PRIVATE 
                                       jabber
                                       cxxopts::cxxopts)
   list(APPEND JABBER_APP_TARGETS jabber_plot)                                       
endif()

# Build Jabber PSD app
if (JABBER_BUILD_PSD)
   add_executable(jabber_psd jabber_psd.cpp)
   target_link_libraries(jabber_psd PRIVATE 
                                       jabber
                                       cxxopts::cxxopts)
   target_include_directories(jabber_psd PRIVATE
                     ${CMAKE_SOURCE_DIR}/externals/pocketfft)
   list(APPEND JABBER_APP_TARGETS jabber_psd)
endif()

install(TARGETS ${JABBER_APP_TARGETS}
   DESTINATION ${CMAKE_INSTALL_BINDIR})
