# This makefile can be used to compile iran3d.

# For Linux or UNIX systems I recommend using the gcc for compilation.
# On Windows I recommend using mingw32.
#CC = gcc
#CC = mingw32-gcc
CC = mpicc

#CFLAGS = -O2 -Wall -ansi -pedantic
#CFLAGS = -O2 -Wall -pedantic
#CFLAGS = -O1 -Wall
CFLAGS = -O1 -Wall -m32# 32bit, -arch i386
LDFLAGS = -lm -m32# -v
#LIBFLAG = -L. -lstruct_s# serial
LIBFLAG = -L. -lstruct_m# mpi

# Notes on warning level:
# Using -Wall and -pedantic will return many warnings, because of non-allowed
# comment styles in the codes. I recommend using just -Wall.

# Notes on the compiler optimization options -OX:
# It is not recommended to use -O2 or -O3. On Linux systems no problems using
# these options have been observed so far; however, on windows systems porgram
# crashes and hang-ups did occur when using -O2 or -O3.
# In windows, the inverse sqr tables will not work with -O2.
# If unsure, compile all with -O1 only.

im3d : im3d.o mpimod.o const.o init.o material.o target.o matrix.o index.o magic.o fileio.o aivxyz.o cfgwriter.o mshwriter.o bulk.o utils.o random.o
	$(CC) $(LDFLAGS) -o im3d im3d.o mpimod.o const.o init.o material.o target.o matrix.o index.o magic.o fileio.o aivxyz.o cfgwriter.o mshwriter.o bulk.o utils.o random.o $(LIBFLAG)

im3d.o: im3d.h im3d.c

mpimod.o: mpimod.h mpimod.c

const.o: const.h const.c

init.o: init.h init.c

material.o: material.h material.c

target.o: target.h target.c

matrix.o: matrix.h matrix.c

index.o: index.h index.c

magic.o: magic.h magic.c

fileio.o: fileio.h fileio.c

aivxyz.o: aivxyz.h aivxyz.c

cfgwriter.o: cfgwriter.h cfgwriter.c

mshwriter.o: mshwriter.h mshwriter.c

bulk.o: bulk.h bulk.c

utils.o: utils.h utils.c

random.o: random.h random.c

.c.o :
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f im3d *.o
clear:
	rm *.o
