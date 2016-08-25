rm snap/*
mpirun -np 8 ~/bin/lmp_mpi-modc  -in lammps.input.sh  #path to your lammps exec

rm t-vs-sf.dat
awk 'BEGIN {print 0,1,0,0,0}' >> t-vs-sf.dat
for i in snap/SF*
do 
	awk '{if(NR==3) {print $2,$6,($6-1),$7,$8}}' "$i" >> t-vs-sf.dat
done 
