# $Id: Makefile,v 1.1 2007/11/23 04:48:52 kichiki Exp $

CC           = cc
#LD           = libtool --mode=link $(CC)
LD           = libtool --mode=link --tag=CC $(CC)

SRC_DIR = ../src
LIBSTOKES_SRC_DIR = ../../libstokes/src

CFLAGS	= -Wall \
	-pipe -O3 -fomit-frame-pointer -ffast-math\
	-march=prescott -mfpmath=sse\
	-I$(SRC_DIR)\
	-I$(LIBSTOKES_SRC_DIR)
# for SUSE-Linux
#CFLAGS	= -Wall \
#	-pipe -O3 -fomit-frame-pointer -ffast-math\
#	-I$(SRC_DIR)\
#	-I$(LIBSTOKES_SRC_DIR)

ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.3 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
	-larpack
# for SUSE-Linux
#ATLAS_LIBS	= \
#	-lgfortran -lblas -llapack\
#	-larpack

LDFLAGS      = 
LIBS         =	../src/libiter.la\
	     	-lm\
		-L/usr/local/lib\
		$(ATLAS_LIBS)

all:		test-libiter

dgetri_c.o:	$(LIBSTOKES_SRC_DIR)/dgetri_c.c
bench.o:	$(LIBSTOKES_SRC_DIR)/bench.c
# for gnu-make
#dgetri_c.o:	$(LIBSTOKES_SRC_DIR)/dgetri_c.c
#	$(CC) $(CFLAGS) -c -o dgetri_c.o $(LIBSTOKES_SRC_DIR)/dgetri_c.c
#bench.o:	$(LIBSTOKES_SRC_DIR)/bench.c
#	$(CC) $(CFLAGS) -c -o bench.o $(LIBSTOKES_SRC_DIR)/bench.c


OBJ =	test-libiter.o\
	check.o	\
	check-iter-gen.o\
	toeplitz.o\
	dgetri_c.o\
	bench.o


test-libiter:	$(OBJ) ../src/libiter.la
	$(LD) $(LDFLAGS) -o test-libiter $(OBJ) $(LIBS)

clean:	
	rm -f *~ *.o test-libiter
