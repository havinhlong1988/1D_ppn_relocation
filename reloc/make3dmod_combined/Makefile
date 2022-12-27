EXE = make3dmod
CC = gcc
FF = gfortran
INCDIR = ../
CFLAGS = -O3 -I$(INCDIR)
FFLAGS = -O3 -mcmodel=medium -I$(INCDIR)
USER_OBJ = make3dmod.o 3dmod.params

.PHONY: run clean
run: $(EXE) 3dmod.params
	./$(EXE)

$(EXE): $(USER_OBJ)
	$(FF) $(FFLAGS) -o $@ *.o

%.o: %.c
	$(CC) $(CFLAGS) -c $<<
%.o: %.f90
	$(FF) $(FFLAGS) -c $<
	
clean:
	rm -f *.o $(EXE) MOD.ini MOD3D

