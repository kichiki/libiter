# Build a libtool library, libiter.la for installation in libdir.
# $Id: Makefile.am,v 1.2 2006/09/26 18:06:07 ichiki Exp $
lib_LTLIBRARIES = libiter.la

libiter_la_SOURCES = \
  myblas.c \
  gmres.c \
  bi-cgstab.c \
  orthomin.c

libiter_la_LDFLAGS = -version-info 0:0:0

include_HEADERS = libiter.h

# CFLAGS = -Wall -O3 -mtune=pentium
AM_CFLAGS = -Wall -O3 -mtune=pentium
