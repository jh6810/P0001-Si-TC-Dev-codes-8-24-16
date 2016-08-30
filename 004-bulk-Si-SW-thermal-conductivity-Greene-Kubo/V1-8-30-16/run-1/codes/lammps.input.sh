##NOTE: THIS IS NOT A SHELL SCRIP IT IS A LAMMPS INPUT FILE THE EXT .sh 
##      IS ONLY USED SO THAT TEXT EDITORS TREAT IT LIKE A SHELL SCIPT 
##      BECAUSE BOTH HAVE THE SAME COMMENT FORMALISM i.e. # 

##--------------------------------------------------------------------
##CONTACT 
##--------------------------------------------------------------------
## JAMES HICKMAN (GMU) 
## 8-23-2016

##--------------------------------------------------------------------
##COMMENTS  
##--------------------------------------------------------------------

##SEE SHELL SCRIPT 'run.sh' FOR RUN COMMAND  


##--------------------------------------------------------------------
##useful variables 
##--------------------------------------------------------------------

##using default metal timestep: 1000 steps=1ps   1,000,000 steps= 1 ns

##0K PARAM 
variable a_o equal 5   	#Guess for un-min lattice constant
variable rep_x equal 5   	#repeat lattice this many times in given direction
variable rep_y equal 10  	#repeat lattice this many times in given direction
variable rep_z equal 5   	#repeat lattice this many times in given direction
variable T equal 500  		#desired temperature of run

##NPT PARAM 
variable N_EQ equal 2500       #number of steps in equilibration period (NVT AND NPT)
variable N_PROD_NPT equal 10000    #number of steps in production period 
variable N_PROD_NVT equal 2500    #number of steps in production period 
variable N_GK equal 1000000  #1000 steps=1ps   1,000,000 steps= 1 ns

variable thermo_save equal 100   #TD sampling rate 

variable LOC1 string ../output/output-1/
variable LOC2 string ../output/output-2/

###--------------------------------------------------------------------
###SIMULATION-1: 0K MINIMIZATION RUN
###--------------------------------------------------------------------

#log  ${LOC1}/structure-gen.log
log temporary.log

units 		metal
dimension	3
boundary	p	p	p
atom_style	atomic

#lammps commands diamond,fcc,hcp all construct in the 
#conventional sense (simple cubic lattice)
lattice diamond ${a_o} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1  
region box block 0 1 0 1 0 1 units lattice
create_box   1 box
create_atoms 1 box
replicate ${rep_x} ${rep_y} ${rep_z}

##SOME QUANTITIES OF INTEREST
variable a equal (lx)/(${rep_x}) #INSTANTAEOUS LATTICE CONSTANT 
compute 1 all pe/atom

compute 2 all reduce ave c_1

###load potential
pair_style sw
pair_coeff * * pot/Si.sw Si
mass * 28.0855

##RELAX STRUCTURE
fix 1 all box/relax iso 0.0 vmax 0.001
thermo 1 
thermo_style custom step pe c_2 press pxx pyy pzz v_a
min_style cg 
minimize 1e-25 1e-25 5000 10000 

write_data lammps_data.0K.minimized #data file

print "minimum-lattice constants: ${a}"

##DEFINE EQ 0K BLOCK SIZE FOR THERMAL EXPANSION CALCULATION
variable lxo equal ${rep_x}*${a}  #treated as constants 
variable lyo equal ${rep_y}*${a}
variable lzo equal ${rep_z}*${a}
clear

####--------------------------------------------------------------------
####SIMULATION-2: Run NPT to determine thermal expansion factor for chosen temperature 
####--------------------------------------------------------------------
#log temporary.log

#units 		metal
#dimension	3
#boundary	p	p	p
#atom_style	atomic

##load structure
#read_data lammps_data.0K.minimized

#velocity all create ${T} 4928459 rot yes mom yes dist gaussian

##SCALE FACTORS
#variable sf equal (lx/${lxo}+ly/${lyo}+lz/${lzo})/3.0

##load potential
#pair_style sw
#pair_coeff * * pot/Si.sw Si
#mass * 28.0855


#reset_timestep 0
#thermo ${thermo_save}
#thermo_style custom time temp press v_sf
#fix 1 all npt temp ${T} ${T} 0.1 aniso 0.0 0.0 1.0 

##EQ RUN 
#run ${N_EQ} 

##PROD RUN 
#reset_timestep 0
#variable Ti equal temp 	#K
#variable Pi equal press*0.0001 #GPa  

##Nevery=2, Nrepeat=6, and Nfreq=100, then values on timesteps 90,92,94,96,98,100
#fix ave_SF all ave/time 1 ${N_PROD_NPT} ${N_PROD_NPT} v_sf 

#run ${N_PROD_NPT} 

#variable ave_SF_temp equal f_ave_SF
#variable ave_SF equal ${ave_SF_temp} # freezes value using ${}, otherwise lammps will try to compute using the fix 
#print ${ave_SF}

##write_data ${LOC1}/NPT-FINAL-STRUCTURE.${T}K #data file

#clear 

#####--------------------------------------------------------------------
#####SIMULATION-3: Scale 0K structure to relevant temperature and
#####              randomize positions in NVT (also check pressure state)
#####--------------------------------------------------------------------

##ASSUMING SCALE FACTOR FOR CORRECT TEMP IS KNOWN YOU CAN SKIP TO HERE
variable ave_SF equal 1.00194017060324

log temporary.log

units 		metal
dimension	3
boundary	p	p	p
atom_style	atomic

##load structure
read_data lammps_data.0K.minimized

##scale block and randomize velocities
print ${ave_SF}

change_box all x scale ${ave_SF} y scale ${ave_SF} z scale ${ave_SF} remap
velocity all create ${T} 4928459 rot yes mom yes dist gaussian

#load potential
pair_style sw
pair_coeff * * pot/Si.sw Si
mass * 28.0855

thermo ${thermo_save}
thermo_style custom time temp press 
fix 1 all nvt temp ${T} ${T} 0.1 

##EQ RUN 
run ${N_EQ} 
#log output-2/log.TC.${T}  #save TD output to this file	

##PROD RUN (just to check the average pressure) 
reset_timestep 0

variable Pxx equal pxx*0.0001 #GPa  
variable Pyy equal pyy*0.0001 #GPa  
variable Pzz equal pzz*0.0001 #GPa  
variable Pxy equal pxy*0.0001 #GPa  
variable Pxz equal pxz*0.0001 #GPa  
variable Pyz equal pyz*0.0001 #GPa  

##Nevery=2, Nrepeat=6, and Nfreq=100, then values on timesteps 90,92,94,96,98,100
fix ave_press all ave/time 1 ${N_PROD_NVT} ${N_PROD_NVT} v_Pxx v_Pyy v_Pzz v_Pxy v_Pxz v_Pyz  file  ${LOC1}/NVT-ave-pressure.dat

run ${N_PROD_NVT} 

write_data ${LOC1}/NVT-EQ-STRUCTURE.lammps #data file

clear 

#####--------------------------------------------------------------------
#####SIMULATION-4: Greene Kubo
#####--------------------------------------------------------------------


variable thermo_save equal 1 
variable snap_save equal ${N_GK}

#variable data_save equal ${N_steps}

units 		metal
dimension	3
boundary	p	p	p
atom_style	atomic

##READ ATOMS
read_data ${LOC1}/NVT-EQ-STRUCTURE.lammps
velocity all create ${T} 4928459 rot yes mom yes dist gaussian

##dump dump1 all custom ${snap_save} snap/snap.*.lammps id type x y z 

#load potential
pair_style sw
pair_coeff * * pot/Si.sw Si
mass * 28.0855

compute      myKE all ke/atom
compute      myPE all pe/atom
compute      myStress all stress/atom NULL virial
compute      flux all heat/flux myKE myPE myStress
variable     Jx equal c_flux[1]/vol
variable     Jy equal c_flux[2]/vol
variable     Jz equal c_flux[3]/vol
variable     J equal (c_flux[1]+c_flux[2]+c_flux[3])/(3*vol)
variable     p_gPa equal 0.0001*press

variable Pxx equal pxx*0.0001 #GPa  
variable Pyy equal pyy*0.0001 #GPa  
variable Pzz equal pzz*0.0001 #GPa  
variable Pxy equal pxy*0.0001 #GPa  
variable Pxz equal pxz*0.0001 #GPa  
variable Pyz equal pyz*0.0001 #GPa  

##Nevery=2, Nrepeat=6, and Nfreq=100, then values on timesteps 90,92,94,96,98,100
fix ave_press all ave/time 1 ${N_GK} ${N_GK} v_Pxx v_Pyy v_Pzz v_Pxy v_Pxz v_Pyz  file  ${LOC2}/NVT-ave-pressure-GK.dat 

variable V equal vol 
print ${V} file ${LOC2}/volume.dat
log ${LOC2}/NVE.log

fix 1 all nve

thermo ${thermo_save}
thermo_style custom time v_J temp v_Jx v_Jy v_Jz press etotal ke pe
run ${N_GK} 

write_data ${LOC2}/data.${N_GK}


