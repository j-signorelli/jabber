set(JABBER_CORE_SOURCES
      interpolant.cpp
      acoustic_field.cpp
      transfer_functions.cpp
      kernels.cpp
      psd.cpp)

set(JABBER_CORE_HEADERS
      interpolant.hpp
      acoustic_field.hpp
      transfer_functions.hpp
      kernels.hpp
      psd.hpp
      core.hpp)

if(JABBER_ENABLE_OPENMP)
   target_link_libraries(jabber PRIVATE OpenMP::OpenMP_CXX)
endif()

append_relative_paths(JABBER_CORE_SOURCES)
append_relative_paths(JABBER_CORE_HEADERS)

set(JABBER_SOURCES ${JABBER_SOURCES} ${JABBER_CORE_SOURCES} PARENT_SCOPE)
set(JABBER_HEADERS ${JABBER_HEADERS} ${JABBER_CORE_HEADERS} PARENT_SCOPE)
