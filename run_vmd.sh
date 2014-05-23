#! /bin/bash
#$ -S /bin/bash
#$ -R yes
#$ -V
#$ -cwd
#$ -N cluster
#$ -j y
#$ -q cpu_long
#$ -P kenprj

#module load gcc/4.7.0

/vault/tools/vmd/bin/vmd  -dispdev text -eofexit <  /cavern/alberto/Ffield/OPENMM/Scripts/analyze_phi.tcl -args Data/trajectory.pdb
