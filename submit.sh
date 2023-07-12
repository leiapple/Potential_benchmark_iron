#!/bin/bash
#------------------------------
# This is a main/controlling code written in SHELL.
# Please read through README.md file carefully before using the script.
# DATE: 1 Dec 2021
# AUTHOR: lei.zhang@rug.nl
#-------------------------------
# 1st updates: 11 Oct 2022: New function
# Major change: all calculations are submitted by SLURM.
#-------------------------------
# 2nd updates: 2 Dec 2022 
# Major change: 
# - separate the lammps calculation and postprocessing
# - copy input sript and perform the calculation in one folder, 
# which allows multiple tests without replicating folders manually.
# - fix the inconsistency of Bain path calculation
# - Merge the surface energy input into one script
#------------------------------
# 3rd updates: 10 June 2023
# - add the calculation of T-S curve
# - add jace1x potential
#------------------------------
#SBATCH --job-name=IAP_test
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=parallelshort
#SBATCH --time=2:00:00
#SBATCH --error=slurm-%j.stderr
#SBATCH --output=slurm-%j.stdout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lei.zhang@rug.nl

# clear caches
rm dump*
rm *.csv
rm sfe*
rm -r ./data
rm -r ./plots
rm in.*
rm *.mod
rm *.py
rm results.txt
rm *.log
rm potential/potential.in
#**********************************
# Customize section
#**********************************
# your email address
eaddress="lei.zhang@rug.nl"
# load the dependent modules of the cluster.
module load OpenMPI/4.1.4-GCC-11.3.0
module load Python/3.10.4-GCCcore-11.3.0
# Full path to LAMMPS excutable. 
LMMP="/home1/p301616/software/lammps_jace/build/lmp"
# the foldername of lammps input scripts (with fullpath)
lmp_inps="/home1/p301616/pot_testing/lmps_inputs"
pps_python="/home1/p301616/pot_testing/py_pps"
# Choose interatomic potential
pot="ace"
# input the potential file name
potfilename=`ls ${PWD}/potential/`

# For GAP, one need to define the right coeff
if [[ ${pot} = "ace" ]]; then
	pstyle="pace"
	pcoeff="${PWD}/potential/${potfilename} Fe"
elif [[ ${pot} = "gap" ]]; then
	pstyle="quip"
	pcoeff="${PWD}/potential/${potfilename} 'IP GAP' 26"
elif [[ ${pot} = "n2p2" ]]; then
	pstyle="hdnnp 6.5 dir ${PWD}/potential showew no showewsum 1000 resetew yes maxew 1000000 cflength 1 cfenergy 1"
	pcoeff="Fe"
elif [[ ${pot} = "ann" ]]; then
	pstyle="aenet"
	pcoeff="v-1 Fe 15tw-15tw.nn Fe"
elif [[ ${pot} = "mtp" ]]; then
	pstyle="mlip ${PWD}/potential/mlip.ini"
	poceff=""
fi

# Generate the interatomic potential file
cat >./potential/potential.in <<EOF
# Define the interatomic potential
pair_style ${pstyle}
pair_coeff * * ${pcoeff} 
EOF
# create a data folder
mkdir data

#**********************************
# Get the information 
#**********************************
## locate the folder and grep the folder name
fullpath=${PWD}
potential_name=`echo $(basename $fullpath)`
# Grep the potential version and echo to results file
echo '#**********************************' | tee -a  ./data/results.txt
echo 'Potential basis set:' ${potential_name} | tee -a ./data/results.txt
awk '/^pair_style*/' ./potential/potential.in | tee -a ./data/results.txt
awk '/^pair_coeff*/' ./potential/potential.in | tee -a ./data/results.txt
echo '#**********************************' | tee -a ./data/results.txt

#**********************************
# Calculation section
#**********************************
# E-V curve 
cp ${lmp_inps}/in.eos .
mpirun -n 4 ${LMMP} -in in.eos -v folder ${potential_name}
# fit EOS
cp ${pps_python}/eos-fit.py .
python eos-fit.py
cp volume.dat ./data/eos_mlip.csv
# Get lattice parameter
a0=$(grep 'a0 =' ./data/results.txt | awk '{print $3}')

# Vacancy formation energy
cp ${lmp_inps}/in.vac .
srun ${LMMP} -in in.vac -v lat ${a0}

# Calculation of elastic constants.--------------------------------
cp ${lmp_inps}/in.elastic .
cp ${lmp_inps}/*.mod .
srun ${LMMP} -in in.elastic -v lat ${a0}

# Calculation of surface energies.---------------------------------
cp ${lmp_inps}/in.surf* .
# (100) plane
srun ${LMMP} -in in.surf1 -v lat ${a0}
# (110) plane
srun ${LMMP} -in in.surf2 -v lat ${a0}
# (111) plane
srun ${LMMP} -in in.surf3 -v lat ${a0}
# (112) plane
srun ${LMMP} -in in.surf4 -v lat ${a0}

# Bain path calculation.------------------------------------------
cp ${lmp_inps}/in.bain_path .
${LMMP} -in in.bain_path -v lat ${a0}
cp bain_path.csv ./data

# Stacking fault energy---------------------------------------------
cp ${lmp_inps}/in.sfe_* .
srun ${LMMP} -in in.sfe_110 -v lat ${a0}
srun ${LMMP} -in in.sfe_112 -v lat ${a0}
cp ./sfe_110.csv ./data
cp ./sfe_112.csv ./data

# Traction-separatio curve------------------------------------------
cp ${lmp_inps}/in.ts_* .
srun ${LMMP} -in in.ts_100 -v lat ${a0}
srun ${LMMP} -in in.ts_110 -v lat ${a0}
cp ./ts_100.csv ./data
cp ./ts_110.csv ./data

#**********************************
# Plotting section
#**********************************
# Execute python script to do the plots.py------------------------------
# Plot E-V curve and Bain path
cp -r /home1/p301616/pot_testing/REF_DATA . 
mkdir plots
cd plots
cp ${pps_python}/eos_bain.py .
cp ${pps_python}/sfe.py .
cp ${pps_python}/ts.py .
python eos_bain.py
python sfe.py
python ts.py
rm *.py

echo "Finish plotting results!"

# delete all lammps inputs
cd ..
rm in.*
rm *.mod

# Send email of the plots as the attached file.
# Mail the results ---------------------------------------------------
mail -s "Basic Properties of iron predicted by IAP"  -a ./data/results.txt -a ./plots/eos_bp.png -a ./plots/sfe.png "${eaddress}" <<EOF
Please check the performance of interatomic potential: ${potential_name}
EOF
echo "Mail the results successful!"
