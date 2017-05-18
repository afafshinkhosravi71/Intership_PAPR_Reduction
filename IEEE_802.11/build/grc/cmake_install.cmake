# Install script for directory: /home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gnuradio/grc/blocks" TYPE FILE FILES
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_moving_average_xx.xml"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_sync_long.xml"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_ether_encap.xml"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_decode_mac.xml"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_sync_short.xml"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_parse_mac.xml"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_chunks_to_symbols.xml"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_mapper.xml"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_mac.xml"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/grc/ieee802_11_frame_equalizer.xml"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
