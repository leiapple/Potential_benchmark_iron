# Nudged elastic band calculations of Peierls barrier in LAMMPS 

This folder gives an example of the how the NEB is carried out in [publication](https://arxiv.org/abs/2307.10072). Before performing NEB calculations, the initial and final replicas of the system need to be prepared. In the case of dislocations, the initial and final image differs by a glide distance of Burgers vector. In this example, we use the screw dislocation in iron to demonstrate the workflow. One can start with any suitable potential for any bcc systems one is interested. Please set the potential parameters in any LAMMPS script before usage.

1. Prepare initial dislocation structure.
In folder `01_Initial_Configuration`, run the lammps script `input_BCC_screw`. This script will generate the equilibrium dislocation structure in a periodic array of dislocation (PAD) configuration. The script adopts a mixed usage of both conjugate gradient and FIRE minimizer to avoid local minima of the potential energy surface. Moreover, an incremental shear force is applied to the dislocation until it moves to the next Peierls valley. 
2. Prepare the final image of the NEB calculation.
Copy the restart file and the moved dislocation dump file to folder `02_Final_Configuration`, change lammps script `input_final_BCC`accordingly. Run the lammps script.
3. Run the bash script `run.sh` in folder `03_NEB_Input_Data` to get the initial and final replica for the NEB calculation in lammps.
4. Run `input_neb` to get the NEB path. 
5. Run python script `plot_peierls.py` to get the Peierl barrier energy vs. reaction coordinate.

Please note that this example is meant for users with LAMMPS experience. Please do not hesitate to contact `lei.zhang@rug.nl` if you encounter any problems. 