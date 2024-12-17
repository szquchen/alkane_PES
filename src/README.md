## bemsa*.f90 files
These contain the permutationally invariant polynomials (PIPs), used as the basis functions in this potential energy surface.

Any change in these files is not recommended.

## example.f90
Sample main program. Provides an example showing how to call the potential. to run the example, download the package and compile the program using the Makefile provided (with proper modifications based on the compiler), and copy the resulting executable (example.x), the coeff.dat, and test.xyz to the same directory. Run\
./example.x

## pes_shell.f90
Interface if the users want to link the potential to their driver program, such as molecular dynamics.

The key subroutines/functions in this file are described below:
(1) pes_init(): reading in the linear coefficients. Must be called once before calling the potential
(2) subroutine getpot(x1d, pot): calculates the potential energy. x1d(3*natm) is the Cartesian coordinate in
    bohr, and the pot is the output, and the energy is in unit hartree
(3) subroutine pot_grad(x1d, pot, grad1d): calculate the potential and gradient simultaneously (which will be slightly faster
    than calculating them separately). x1d and pot are the same as those in getpot(), and grad1d is the calculated gradient,
    in hartree/bohr
