#!/bin/bash

# The script is used to submit k-test jobs for bcc metals.

# Full path to LAMMPS excutable. 
LMMP="/home1/p301616/software/lammps/build_gap/lmp"
# the foldername of lammps input scripts (with fullpath)
lmp_inps="/home1/p301616/pot_testing/lmps_inputs"

# first, solve the coeffs. 
module load MATLAB
matlab -nodisplay -r Solve_aniso_coeff_v2

# get the equilibrium constants
a0=$(grep 'a0 =' ./bench/data/results.txt | awk '{print $3}')
mass=95.95

for crksys in $(seq 1 1 4)
do

KI=$(grep '#K_I=' ./lefm_coeffs/lefm_paras.CrackSystem_${crksys} | awk '{print $2}')
Kstart=`printf "%.0f" $(bc <<< "$KI*100-10")`
Kstop=`printf "%.0f" $(bc <<< "$Kstart+100")`

mkdir CrackSystem_${crksys}
cd CrackSystem_${crksys}
cp ${lmp_inps}/in.cracksystem_${crksys} .

cat << EOF > submit.sh
#!/bin/bash
#SBATCH --job-name=crksys.${crksys}
#SBATCH --partition=regularmedium
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --error=slurm-%j.stderr
#SBATCH --output=slurm-%j.stdout

module restore gap

srun ${LMMP} -in in.cracksystem_${crksys} -v a0 ${a0} -v m ${mass} -v CrkSys ${crksys} -v Kstart ${Kstart} -v Kstop ${Kstop}
EOF

sbatch submit.sh
cd ..
done 
