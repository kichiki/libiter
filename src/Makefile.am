# Build a libtool library, libiter.la for installation in libdir.
# $Id: Makefile.am,v 1.4 2006/10/05 21:54:31 ichiki Exp $
lib_LTLIBRARIES = libiter.la

libiter_la_SOURCES = \
  libiter.c \
  myblas.c \
  gmres.c \
  bi-cgstab.c \
  orthomin.c

libiter_la_LDFLAGS = -version-info 1:0:0

include_HEADERS = libiter.h

AM_CFLAGS = -Wall -O3 -mtune=pentium
