EXE = make_td_crust
SRCDIR = ../../src/reloc/
CC = gcc
FF = gfortran
CFLAGS = -Wall
FFLAGS =
USER_OBJ = distazsub.o make_td_crust.o func.o 
#USER_DIR = makeinput

.PHONY: run clean backup
run: $(EXE)
#	./$(EXE) < input

$(EXE): $(USER_OBJ)
	$(FF) -o $@ *.o

%.o: %.c para.inp
	$(CC) $(CFLAGS) -c $<
%.o: %.f90 para.inp
	$(FF) $(FFLAGS) -c $<
	
clean:
	rm -f *.o
backup:
	cp -f *.c *.f90 *.inp input readme Makefile $(SRCDIR)
	cp -rf $(USER_DIR) $(SRCDIR)

