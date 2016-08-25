
#--------------------------------------------------------------------
#CONTACT 
#--------------------------------------------------------------------
# JAMES HICKMAN (GMU) 
# 8-23-2016

#--------------------------------------------------------------------
#DESCRIPTION 
#--------------------------------------------------------------------

#SEE SHELL SCRIPT 'run.sh' FOR RUN COMMAND  

#LAMMPS CODE TO COMPUTE THERMAL EXPANSION CURVES AS A FUNCTION OF
#OF TEMPERATURE IN SINGLE COMPONENT SYSTEM

#NOTE-1: CURRENT FORM IS WRITEN FOR Cu BUT HAS BEEN TESTED WITH Si AS WELL

#NOTE-2: SHOULD BE ABLE TO USE FOR ANY POTENTIAL WITH MINIMIAL 
#	 MODIFICATIONS TO THE CODE BELOW

#--------------------------------------------------------------------
#useful variables 
#--------------------------------------------------------------------

#using default metal timestep: 1000 steps=1ps   1,000,000 steps= 1 ns

#0K PARAM 
variable rep equal 5		#NUM of replications of block in cubic directions
variable a_o equal 3.6   	#Guess for un-min lattice constant

#NPT PARAM 
variable N_EQ equal 5000        #number of steps in equilibration period 
variable N_PROD equal 100000    #number of steps in production period 
variable thermo_save equal 25   #TD sampling rate 
variable DT equal 100	 	#temperature incriment in K 	
variable NUM_TEMPS equal 14	#temp mesh 1*DT,2*DT....NUM_TEMPS*DT	

#--------------------------------------------------------------------
#0K MINIMIZATION RUN
#--------------------------------------------------------------------

units 		metal
dimension	3
boundary	p	p	p
atom_style	atomic

#lammps commands diamond,fcc,hcp all construct in the 
#conventional sense (simple cubic lattice)
lattice fcc ${a_o} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1  
region box block 0 1 0 1 0 1 units lattice
create_box   1 box
create_atoms 1 box
replicate ${rep} ${rep} ${rep}

#SOME QUANTITIES OF INTEREST
variable a equal (lx+ly+lz)/(3.0*${rep}) #INSTANTAEOUS LATTICE CONSTANT 
compute 1 all pe/atom
compute 2 all reduce ave c_1

##load potential
pair_style	eam/alloy
pair_coeff	* * pot/cu_mishin_2001.eam.alloy  Cu


#RELAX STRUCTURE
fix 1 all box/relax iso 0.0 vmax 0.001
thermo 1 
thermo_style custom step pe c_2 press pxx pyy pzz v_a
min_style cg 
minimize 1e-25 1e-25 5000 10000 

write_data snap/lammps_data.0K.final #data file

print "minimum-lattice constants: ${a}"

#DEFINE EQ 0K BLOCK SIZE FOR THERMAL EXPANSION CALCULATION
variable lxo equal ${rep}*${a}
variable lyo equal ${rep}*${a}
variable lzo equal ${rep}*${a}
print "Lx: ${lxo}"
clear

#--------------------------------------------------------------------
#NPT-loop for averaging
#--------------------------------------------------------------------
variable i loop ${NUM_TEMPS} #INCREMENTS i++
label loop #START OF LOOP 
	print "${i} ${lxo}"
	variable T equal ${i}*${DT}
#	#-----------------NVE RUN--------------------
	log snap/eq.log.${T}  #save TD output to this file	

	units 		metal
	dimension	3
	boundary	p	p	p
	atom_style	atomic

	variable Tm1 equal ${T}-${DT} #Tm1=T minus 1

#	#load structure
	if "${i}==1" then "read_data snap/lammps_data.0K.final" else "read_data snap/data.${Tm1}"

	velocity all create ${T} 4928459 rot yes mom yes dist gaussian

	#SCALE FACTORS
	variable sfx equal lx/${lxo}
	variable sfy equal ly/${lyo}
	variable sfz equal lz/${lzo}
	variable sf equal (lx/${lxo}+ly/${lyo}+lz/${lzo})/3.0

	##load potential
	pair_style	eam/alloy
	pair_coeff	* * pot/cu_mishin_2001.eam.alloy  Cu

	reset_timestep 0
	thermo ${thermo_save}
	thermo_style custom time temp press v_sf
	fix 1 all npt temp ${T} ${T} 0.1 aniso 0.0 0.0 1.0 

	#EQ RUN 
	run ${N_EQ} 

	#PROD RUN 
	reset_timestep 0
	log snap/prod.log.${T}  #save TD output to this file	
	variable Ti equal temp 	#K
	variable Pi equal press*0.0001 #GPa  
	#Nevery=2, Nrepeat=6, and Nfreq=100, then values on timesteps 90,92,94,96,98,100
	fix ave1 all ave/time 1 ${N_PROD} ${N_PROD} v_T v_sfx v_sfy v_sfz v_sf v_Ti v_Pi file snap/SF.${T}.log 
	run ${N_PROD} 

	write_data snap/data.${T}

	clear 

	next i

jump SELF loop #JUMP TO POINT LABEL-LOOP IN THIS FILE (SELF) 


