# Targeted Transferable PES for n-alkanes
This is a targeted transferable potential energy surface (PES) for n-alkanes. It is based on the many-body permutationally invariant polynomial approach desecribed in Ref. [1], and is trained on ~250,000 C<sub>14</sub>H<sub>30</sub> energies computed at B3LYP/cc-pVDZ level of theory. It can be applied to linear alkanes, and its fidelity has been tested on alkanes from C<sub>4</sub>H<sub>10</sub> to C<sub>30</sub>H<sub>62</sub>.

# Usage
Source code is the the folder "src", and the parameters of the potential is in coeff.dat.

This is just the potential energy surface; users need to link it to their own driver program, such as a molecular dynamics package. The src/pes_shell.f90 provides the interface for this purpose. We do provide a simple main program "example.f90", however, to show how to call our potential. Detailed description can be found in the README file inside src folder.  

The PES is written in Fortran, so the users need a Fortran compiler (such as gfortran or Intel ifort).

# References
[1] Qu, Chen; Houston, Paul L.; Allison, Thomas; Schneider, Barry I.; Bowman, Joel M. DFT-Based Permutationally Invariant Polynomial Potentials Capture the Twists and Turns of C<sub>14</sub>H<sub>30</sub>. J. Chem. Theory Comput. 2024, 20, 9339–9353
