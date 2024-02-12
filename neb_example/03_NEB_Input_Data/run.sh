#!/bin/bash

# copy the restart file
cp ../01_Initial_Configuration/restart.screw .
cp ../02_Final_Configuration/restart.screw.moved . 

module restore set-gap

~/software/lammps_pace/build/lmp  -in input_restart2data_start >log.initial
echo 'Initial replica generated!'

~/software/lammps_pace/build/lmp  -in input_restart2data >log.final
echo 'Final replica generated!'

python data2neb_input.py final.data neb.final

cp initial.data ../04_NEB_Calculation_gamma_64rep/
cp neb.final ../04_NEB_Calculation_gamma_64rep/
echo 'Copy configs to `04_NEB_Calculation`'