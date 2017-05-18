# Install script for directory: /home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11

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

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "ieee802_11_devel")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ieee802-11/gnuradio/ieee802_11" TYPE FILE FILES
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/build/include/ieee802-11/moving_average_ff.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/build/include/ieee802-11/moving_average_cc.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/api.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/chunks_to_symbols.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/constellations.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/decode_mac.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/ether_encap.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/frame_equalizer.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/mac.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/mapper.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/parse_mac.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/signal_field.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/sync_long.h"
    "/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/include/ieee802-11/sync_short.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "ieee802_11_devel")

