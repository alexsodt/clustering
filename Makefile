CC = g++
F77 = gfortran 
CC2 = gcc
LDFLAGS = -llapack -lblas -lgfortran -lstdc++

FLAGS = -O2 
 

.C.o: 
	$(CC) $(INCLUDE) $(GSLINC) -o $*.o -c $(FLAGS) $*.C

.c.o:
	$(CC2)	$(QHULL) $(COPT) $(CDEFN)   -g	-c	$*.c

.f.o:
	$(F77) -O3 -c $*.f -o $*.o

EXECS = extractNucleotides cluster

all: $(EXECS)


clean: 
	rm  *.o $(EXECS)

extractNucleotides: extractNucleotides.o util.o dcd.o pdb.o alignSet.o util.o mutil.o
	$(CC) -o extractNucleotides extractNucleotides.o dcd.o pdb.o alignSet.o util.o mutil.o $(LDFLAGS) 

cluster: cluster.o medoidLR.o util.o dcd.o alignSet.o pdb.o geometry.o comparison.o
	$(CC) -o cluster cluster.o medoidLR.o util.o dcd.o alignSet.o pdb.o geometry.o comparison.o $(LDFLAGS)
