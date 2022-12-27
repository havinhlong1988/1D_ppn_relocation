EXE = ../s3dloc
CC = gcc
FF = gfortran
INCDIR = ../
CFLAGS = -O3 -I$(INCDIR)
FFLAGS = -O3 -mcmodel=medium -I$(INCDIR)
USER_OBJ = distazsub.o surfpath.o s3dloc.o func.o pbr.o pn_path.o subr.o ../para.inp ../input reloc2inv1d.o

.PHONY: run clean
#run: $(EXE) ../para.inp ../input
#	./$(EXE) < input

$(EXE): $(USER_OBJ)
	$(FF) $(FFLAGS) -o $@ *.o

%.o: %.c
	$(CC) $(CFLAGS) -c $<
%.o: %.f
	$(FF) $(FFLAGS) -c $<
%.o: %.f90
	$(FF) $(FFLAGS) -c $<
	
clean:
	rm -f *.o $(EXE)

