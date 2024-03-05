# Interatomic potential test via LAMMPS

**DATE**: 9 Mar 2023

**AUTHOR**: <lei.zhang@rug.nl>

This repository contains the lammps script and python script to calculate the following properties for BCC iron. 
It can be generalised to test other materials as well. 
The purpose is to use it on high performance cluster in case that local machine cannot handle the heavy computation, e.g., machine learning potential. 
With minor modifications, the current workflow can also be used on local machines. 

* Energy-Volume curve
* Elastic constants
* Vacancy formation energy
* Surface energies
* Bain path
* Stacking fault (curves) (110)\[1-11\]\[-112\]
* T-S curve
* Dislocation properties (note that this feature does not integrated in the workflow yet)

## Prerequisites 

1. Installation of LAMMPS with the relevant packages compiled.
2. Installation of PYTHON package: matplotlib, pandas, numpy, scipy.

## How it works?

The workflow is governed by the bash file (submit.sh) which is a slurm job file.
This file copy the lammps inputs and run it, which means that it depends on the system configuration.
Thus, one need to specify:

* email address: eaddress
* load the dependent modules (on cluster): 
* give the lammps excutable: LMMP
* specify the path to lammps and python inputs: lmp_inps & pps_python (Those inputs are provided from the repository).
* The DFT reference data if the one wants to generate the plots (the current data is for BCC iron)
* choose the potential: (five potentials are currently available. ACE, GAP, N2P2, ANN, MTP). To customize your own potential, one needs to define it in the script.

After the proper configuration, one is able to perform potential test.

1. copy folder 'pot_testing' and bash file 'submit.sh' to the cluster.
2. Customize the SLURM job submission parameter.
3. make sure the configuration is corrected.
4. run the test by submiting the job script `sbatch submit.sh`

## For testing the dislocation properties.

Currently the workflow does not include the dislocation benchmark. Everything regarding dislocation calculations are put in the `dislocation` folder, including
* Script to create dislocation using [atomsk](https://atomsk.univ-lille.fr/)
* LAMMPS script to crete various dislocations
For more detials, please refer to another readme file inside dislocaiton folder.
