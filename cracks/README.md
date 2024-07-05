# Fracture K-test 

This test performs lammps calculation for four crack systems, including (100)[010], (100)[001], (110)[001]. (110)[1-10].

1. To perform the K test, a lammps source code should be updated first. One can replace the displace_atoms.cpp under `src` folder in orginal lammps code with the one provided here (`lammps_src`). Then compile lammps as usual.

2. First run the basic test of the potential to get the elastic constant and surface energy. With the basic test, a matlab script is provided to solve the linear elastic fracture mechanics problem (Solve_aniso_coeff_v2.m).

3. Then one is able to run the K test, the lammps inputs are given in `pot_testing/lmps_inputs`

For any further questions, please contact <lei.zhang@rug.nl>.