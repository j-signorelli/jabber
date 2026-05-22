set(JABBER_APP_SOURCES
      config_input.cpp
      common.cpp)

set(JABBER_APP_HEADERS
      params.hpp
      config_input.hpp
      common.hpp
      app.hpp)

target_link_libraries(jabber PRIVATE $<BUILD_INTERFACE:toml11::toml11>)

append_relative_paths(JABBER_APP_SOURCES)
append_relative_paths(JABBER_APP_HEADERS)

set(JABBER_SOURCES ${JABBER_SOURCES} ${JABBER_APP_SOURCES} PARENT_SCOPE)
set(JABBER_HEADERS ${JABBER_HEADERS} ${JABBER_APP_HEADERS} PARENT_SCOPE)
