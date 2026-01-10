INSTALLATION INSTRUCTIONS

This document explains how to compile and run the molecular dynamics program provided in this repository.

Requirements:

              Programming language: Fortran 90
              Compiler: gfortran

Compilation: 

              Open a terminal and navigate to the final_work directory (where verlet.f90 is located).
              Compile the program using gfortran. You can choose any name for the executable. For example:

              gfortran -o md_program verlet.f90



Input Files:

              The program requires two input files, which must be placed in the same directory as the executable:
              1) dat.inp – contains the number of atoms and their initial coordinates and masses.
              2) pot.inp – contains the Lennard-Jones potential parameters (eps and sigma).



Running the Program:
                    
                     Run the program from the terminal:   ./md_program

                     The program will produce the following output files:

                                                                          - coord.out – initial coordinates of atoms.
                                                                          - distances.out – distances between each pair of atoms.
                                                                          - trajectory.xyz – trajectory of the system for all time steps.
                                                                          - Screen output includes total potential, kinetic, and total energies.



Visualization:

              The trajectory.xyz file can be opened using Molden: molden trajectory.xyz
              Adjust visualization settings and animation speed as desired.




NOTES:  

      All modules are included in verlet.f90

