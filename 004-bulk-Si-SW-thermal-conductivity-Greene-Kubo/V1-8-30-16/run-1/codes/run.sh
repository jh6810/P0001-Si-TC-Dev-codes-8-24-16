rm -rf ../output 
mkdir ../output
mkdir ../output/output-1
mkdir ../output/output-2

mpirun -np 8 ~/bin/lmp_mpi-30jul16-misc  -in lammps.input.sh  #path to your lammps exec
