# Build a libtool library, libiter.la for installation in libdir.
# $Id: Makefile.am,v 1.8 2007/11/25 18:47:13 kichiki Exp $
lib_LTLIBRARIES = libiter.la

libiter_la_SOURCES =	\
	libiter.c	\
	myblas.c	\
	gmres.c		\
	bi-cgstab.c	\
	orthomin.c	\
	steepest.c	\
	cg.c		\
	bicg.c		\
	bico.c		\
	qmr.c		\
	cgs.c		\
	bicgstab.c	\
	atpres.c	\
	cgne.c		\
	cg_.c

libiter_la_LDFLAGS = -version-info 3:0:0

include_HEADERS = libiter.h
