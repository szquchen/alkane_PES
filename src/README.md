## bemsa*.f90 files
These contain the permutationally invariant polynomials (PIPs), used as the basis functions in this potential energy surface.

Any change in these files is not recommended.

## example.f90
Sample main program. Provides an example showing how to call the potential. to run the example, download the package and compile the program using the Makefile provided (with proper modifications based on the compiler), and copy the resulting executable (example.x), the coeff.dat, and test.xyz to the same directory. Run\
./example.x

## pes_shell.f90
Interface for linking the potential to the driver program, such as a molecular dynamics package.

The key subroutines in this file are described below:\
(1) pes_init(): reading in the linear coefficients. Must be called once before calling the potential.\
(2) subroutine getpot(x1d, pot): calculates the potential energy, for applications where only the potential energy is needed. The input x1d(3*Natm) is the Cartesian coordinate in bohr; atoms must be in the order C C ... C H H ... H. The output is the energy in hartree.\
(3) subroutine pot_grad(x1d, pot, grad1d): calculate the potential and gradient simultaneously, for applications requiring gradients or forces, such as molecular dynamics. The x1d and pot are the same as those in getpot(), and grad1d is the calculated gradient, in hartree/bohr.\
(4) Other helper subroutines include f_switch(): a switching function to smoothly turn the potential to zero when internuclear distance is larger than the cutoff; hessian(): computes the Hessian matrix using finite difference of the gradient, for applications that require Hessian such as normal-mode analysis, geometry optimization with Newton's Method, etc.
