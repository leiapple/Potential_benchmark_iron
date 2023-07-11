# NOTE: This script can be modified for different atomic structures, 
# units, etc. See in.elastic for more info.
#

# Define the finite deformation size. Try several values of this
# variable to verify that results do not depend on it.
variable up equal 1.0e-5
 
# Define the amount of random jiggle for atoms
# This prevents atoms from staying on saddle points
variable atomjiggle equal 1.0e-5

# Uncomment one of these blocks, depending on what units
# you are using in LAMMPS and for output

# metal units, elastic constants in eV/A^3
#units		metal
#variable cfac equal 6.2414e-7
#variable cunits string eV/A^3

# metal units, elastic constants in GPa
units		metal
variable cfac equal 1.0e-4
variable cunits string GPa

# Define minimization parameters
variable etol equal 0.0 
variable ftol equal 1.0e-10
variable maxiter equal 100
variable maxeval equal 1000
variable dmax equal 1.0e-2

# generate the box and atom positions using a diamond lattice
variable latparam equal ${lat}
variable        box_length equal 2 # ------half box size in lattice parameter units

variable        xdim_1 equal -1*(${latparam}*${box_length})*sqrt(1)-0.001
variable        ydim_1 equal -1*(${latparam}*${box_length})*sqrt(1)-0.001
variable        zdim_1 equal -1*(${latparam}*${box_length})*sqrt(1)-0.001
variable        xdim_2 equal (${latparam}*${box_length})*sqrt(1)
variable        ydim_2 equal (${latparam}*${box_length})*sqrt(1)
variable        zdim_2 equal (${latparam}*${box_length})*sqrt(1)
boundary		p p p
log             elastic.log
lattice         bcc ${latparam} orient x 0 0 1 orient y 1 0 0 orient z 0 1 0
region			box block ${xdim_1} ${xdim_2} ${ydim_1} ${ydim_2} ${zdim_1} ${zdim_2} units box
create_box		1 box
create_atoms	1 box

change_box 		all triclinic
# Need to set mass to something, just to satisfy LAMMPS
mass 			1 1.0e-20
