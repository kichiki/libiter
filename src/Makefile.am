# Build a libtool library, libiter.la for installation in libdir.
# $Id: Makefile.am,v 1.3 2006/09/28 04:28:32 kichiki Exp $
lib_LTLIBRARIES = libiter.la

libiter_la_SOURCES = \
  libiter.c \
  myblas.c \
  gmres.c \
  bi-cgstab.c \
  orthomin.c

libiter_la_LDFLAGS = -version-info 0:0:0

include_HEADERS = libiter.h

AM_CFLAGS = -Wall -O3 -mtune=pentium
