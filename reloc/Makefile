EXE = reloc
SRCDIR = ../../src/reloc/
CC = gcc
FF = gfortran
CFLAGS = -Wall
FFLAGS =
USER_OBJ = func.o subr.o pbr.o pn_path.o s3dloc.o distazsub.o surfpath.o para.inp input reloc2inv1d.o
USER_DIR = makeinput

.PHONY: run clean backup
run: $(EXE) input
	./$(EXE) < input

$(EXE): $(USER_OBJ)
	$(FF) -o $@ *.o

%.o: %.c para.inp
	$(CC) $(CFLAGS) -c $<
%.o: %.f90 para.inp
	$(FF) $(FFLAGS) -c $<
	
clean:
	rm -f *.o ray_p* fort* *.err *.cor *.chk *.out s3*.rms s3*.tt *.tt residual_* *_data err_* *.log temp* tmp* 
backup:
	cp -f *.c *.f90 *.inp input readme Makefile $(SRCDIR)
	cp -rf $(USER_DIR) $(SRCDIR)

