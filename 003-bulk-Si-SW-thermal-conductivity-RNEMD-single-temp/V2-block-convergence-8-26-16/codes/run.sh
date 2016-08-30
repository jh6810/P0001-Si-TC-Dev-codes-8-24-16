LOC1=${PWD} 
cd ../
rm -rf output
mkdir output
cd output
LOC2=${PWD} 
mkdir output-1  #structure gen 
mkdir output-2 	#Thermal conductivity calculations
mkdir plots 	#Thermal conductivity calculations
cd $LOC1




for repxz in 5   
do
for repy in  5 10 15 20 25

do
	cd $LOC1

	num_bin=$(awk -v repy=$repy 'BEGIN {print 2*int(5.4309498*repy/3)}')

	num_atoms=$(awk -v repy=$repy -v repxz=$repxz 'BEGIN {print int(repxz*repxz*repy*8)}')

	num_swap=$(awk -v repy=$repy -v repxz=$repxz 'BEGIN {print int(repxz*repxz*repy*8*100/(1000))}')

	sed "s/REP_C_IN/$repxz/g" lammps.template.sh  > tmp1
	sed "s/NUM_LAYERS_IN/$num_bin/g" tmp1  > tmp2
	sed "s/NUM_ATOMS_IN/$num_atoms/g" tmp2  > tmp3
	sed "s/NUM_SWAP/$num_swap/g" tmp3  > tmp4
	sed "s/REP_Y_IN/$repy/g" tmp4  > lammps.input.sh
	rm tmp1 tmp2 tmp3 tmp4

	mpirun -np 16 ~/bin/lmp_mpi-30jul16-misc  -in lammps.input.sh  #path to your lammps exec
	 

	#POST PROCESS
	cd ${LOC2}/output-2/

	y_mid=$(awk -v rep_y=$repy 'BEGIN {print 5.4309498*rep_y/2}') 
	y_1=$(awk -v rep_y=$repy 'BEGIN {print (5.4309498*rep_y/2)/6}') 
	y_2=$(awk -v rep_y=$repy 'BEGIN {print 5*(5.4309498*rep_y/2)/6}') 

	rm y-vs-T.dat
	rm slope.dat

#	awk -v y_1=${y_1} -v y_2=${y_2} '{if(NR>4) {if($2<(y_2))if($2>y_1){{print $2,$4}}}}' "TC-TEMP-PROFILE-$repxz-$repy" >> y-vs-T.dat

	awk -v y_mid=${y_mid} '{if(NR>4) {if($2<y_mid){print $2,$4}}}' "TC-TEMP-PROFILE-$repxz-$repy" >> y-vs-T.dat

	cp TC-TEMP-PROFILE-$repxz-$repy temp.dat


	cp $LOC1/gnuplot.gnu ./
	gnuplot gnuplot.gnu
	rm temp.dat

	mv plot.png ../plots/$repxz-$repy.png
	dT_dx=$(awk '{print $0}' 'slope.dat') # K/ang
	J=$(awk '{print $0/2.0}' "flux-$repxz-$repy.dat")  # eV/(ps ang^2)

	##J=-k dT/dx
	##1 electron volt / picosecond =1.60217662*10^-7 watts
	##1 ang = 10^-10 m 
	k=$(awk -v J=$J -v dT_dx=$dT_dx 'BEGIN {print (J/dT_dx)*(1.60217*10^(-7))/(10^(-10))}')
	##   eV/(ps*K*ang) to W/(m*K)

	Lx=$(awk -v repxz=$repxz 'BEGIN {print 5.4309498*repxz}')
	Ly=$(awk -v repy=$repy 'BEGIN {print 5.4309498*repy}')

	echo $Lx $Ly $Lx $k >> ../Lx-Ly-Lz-k.dat


done	
done 


