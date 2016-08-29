
# rm -rf output-1 
# mkdir output-1 

# rm -rf output-2 
# mkdir output-2 

# mpirun -np 16 ~/bin/lmp_mpi-30jul16-misc  -in lammps.input.sh  #path to your lammps exec

rm y-vs-T.dat
rm slope.dat
awk '{if(NR>4) {if($2<107){print $2,$4}}}' "output-2/tmp2.profile" >> y-vs-T.dat
gnuplot gnuplot.gnu
dT_dx=$(awk '{print $0}' 'slope.dat') # K/ang
J=$(awk '{print $0/2.0}' 'flux.dat')  # eV/(ps ang^2)

#J=-k dT/dx
#1 electron volt / picosecond =1.60217662*10^-7 watts
#1 ang = 10^-10 m 
k=$(awk -v J=$J dT_dx=$dT_dx'{print (J/dT_dx)*(1.60217*10^(-7))/(10^(-10))}' 'flux.dat')
#   eV/(ps*K*ang) to W/(m*K)
echo $k
