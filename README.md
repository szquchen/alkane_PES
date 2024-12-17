# Targeted Transferable PES for n-alkanes
This is a targeted transferable potential energy surface (PES) for n-alkanes. It is based on the many-body permutationally invariant polynomial approach desecribed in Ref. [1], and is trained on ~250,000 C<sub>14</sub>H<sub>30</sub> energies computed at B3LYP/cc-pVDZ level of theory. It can be applied to linear alkanes, and its fidelity has been tested on alkanes from C<sub>4</sub>H<sub>10</sub> to C<sub>30</sub>H<sub>62</sub>.

# Sample program
The PES is written in Fortran, so the users need a Fortran compiler (such as gfortran or Intel ifort). The source code is in folder "src". The parameters of the PES are written in coeff.dat. Detailed description of the program can be found in the README file inside src folder.

A sample main program "example.f90" is given for illustrating and testing purpose. Compile the program using the Makefile provided (with proper modifications based on the compiler), and copy the resulting executable (example.x), the coeff.dat, and test.xyz to the same directory. Run\
./example.x

Users can also link the PES to their own main driver program by using the functions provided in pes_shell.f90.

# References
[1] Qu, Chen; Houston, Paul L.; Allison, Thomas; Schneider, Barry I.; Bowman, Joel M. DFT-Based Permutationally Invariant Polynomial Potentials Capture the Twists and Turns of C<sub>14</sub>H<sub>30</sub>. J. Chem. Theory Comput. 2024, 20, 9339–9353
