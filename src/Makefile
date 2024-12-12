#compiles pes_fit program
FC = ifort
FFLAGS  = -O3 -r8

OBJ := bemsa2.o  bemsa3.o  bemsa21.o \
       bemsa4.o  bemsa31.o  bemsa22.o \
       pes_shell.o  example.o

example.x : $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o : %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -rf *.o *.mod
