# P0001-LAMMPS-thermal-expansion-curves-8-24-16
# REQUIREMENT: lammps needs to be installed on your computer if this is the case this code can be used 
# to run molecular dynamics to generate thermal expansion curves as a function of temperature for your given potential 

# This example is for an EAM copper potential but the code can be modified very easily for other potentials (it was also tested with a Si Terrsoff potential) 
# to run just modify run.sh to point to your LAMMPS executable and make run.sh  executable chmod a+x run.sh 
# once you execute run.sh the entire process should generate a set of simulations at various temperatures and output a thermal expansion file t-vs-sf.dat
# runtime ~ a few hours using 8 cores on a simple workstation
