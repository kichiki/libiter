# $Id: Makefile,v 1.3 2008/10/18 04:25:22 kichiki Exp $

CC = cc
# for robson gcc version 3
#CC = /usr/bin/gcc
# for robson gcc version 3 (64bit)
#CC = /usr/bin/gcc -m64
# for robson gcc version 4
#CC = /usr/local/bin/gcc4
# for robson gcc version 4 (64bit)
#CC = /usr/local/bin/gcc4 -m64


LD = libtool --mode=link --tag=CC $(CC)
#LD = $(CC)


RYUON_LIB = -L/usr/local/lib
# for robson
#RYUON_LIB = -L/home/ichiki/RYUON3-32/lib
#RYUON_LIB = -L/home/ichiki/RYUON3/lib
#RYUON_LIB = -L/home/ichiki/RYUON4-32/lib
#RYUON_LIB = -L/home/ichiki/RYUON4/lib


LIBITER = ../src/libiter.la
#LIBITER = $(RYUON_LIB) -liter


#OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math
# for FreeBSD
OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math\
	-march=prescott -mfpmath=sse
# for robson
#OPTS = -O3 -mcpu=970 -mtune=970 -maltivec


# for FreeBSD
ATLAS_LIBS = \
	-L/usr/local/lib/gcc-4.2.5 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
	-larpack
#ATLAS_LIBS = \
#	-lg2c -llapack\
#	-larpack
# for SUSE-Linux
#ATLAS_LIBS = \
#	-lgfortran -lblas -llapack\
#	-larpack
# for robson (gcc/g77, 64bit)
#ATLAS_LIBS = -lg2c -L/usr/local/lib64 -lcblas -latlas -L/usr/local/lib64 -larpack
# for robson (gcc/g77, 64bit or 32bit)
#ATLAS_LIBS = \
#	-lg2c \
#	$(RYUON_LIB) -llapack -lblas -larpack

# for robson (gcc4/gfortran, 64bit)
#ATLAS_LIBS = -L/usr/local/lib64 -larpack -lessl -lxlf90_r
# for robson (gcc4/gfortran, 64bit or 32bit)
#ATLAS_LIBS = \
#	-L/usr/local/gcc-4.2.2/lib/gcc/powerpc64-suse-linux/4.2.2 -lgfortran \
#	$(RYUON_LIB) -larpack -llapack -lblas


SRC_DIR = ../src
LIBSTOKES_SRC_DIR = ../../libstokes/src

CFLAGS	= -Wall \
	$(OPTS)\
	-I$(SRC_DIR)\
	-I$(LIBSTOKES_SRC_DIR)

LDFLAGS      = 
LIBS         =	$(LIBITER)\
	     	-lm\
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
	check-cg-pc.o\
	check-gmres-m-pc.o\
	check-gpb-pc.o\
	check-otmk-pc.o\
	check-sta-pc.o\
	check-sta2-pc.o\
	ilu.o\
	bench.o


test-libiter:	$(OBJ)
	$(LD) $(LDFLAGS) -o test-libiter $(OBJ) $(LIBS)

clean:	
	rm -f *~ *.o test-libiter
