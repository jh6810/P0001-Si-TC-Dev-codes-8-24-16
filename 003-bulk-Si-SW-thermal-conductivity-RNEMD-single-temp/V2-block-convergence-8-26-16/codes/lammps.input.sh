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

##MANY STEPS (SIMULATION-1 TRHOUGH 3) COULD BE SKIPPED IF THE THERMAL EXPANSION FACTOR WAS  
##ALREADY KNOWN
 
##--------------------------------------------------------------------
##useful variables 
##--------------------------------------------------------------------

##using default metal timestep: 1000 steps=1ps   1,000,000 steps= 1 ns

##0K PARAM 
variable a_o equal 5   	#Guess for un-min lattice constant
variable rep_x equal 5   	#repeat lattice this many times in given direction
variable rep_y equal 5  	#repeat lattice this many times in given direction
variable rep_z equal 5   	#repeat lattice this many times in given direction
variable T equal 500  		#desired temperature of run

##NPT PARAM 
variable N_EQ equal 7500       #number of steps in equilibration period (NVT AND NPT)
variable N_PROD_NPT equal 10000    #number of steps in production period 
variable N_PROD_NVT equal 2000    #number of steps in production period 
variable thermo_save equal 100   #TD sampling rate 

variable N_TC_EQ equal 50000
variable N_TC equal 200000

#variable N_TC_EQ equal 10000
#variable N_TC equal 100000

variable swap_every equal 100
variable TC_num_layers equal 18 
variable num_swap equal 2 

variable LOC1 string ../output/output-1/
variable LOC2 string ../output/output-2/

###--------------------------------------------------------------------
###SIMULATION-1: 0K MINIMIZATION RUN
###--------------------------------------------------------------------

##log  ${LOC1}/structure-gen.log
#log temporary.log

#units 		metal
#dimension	3
#boundary	p	p	p
#atom_style	atomic

##NOTE: MAKING A VARIABLE CONSTANT IN TIME 
##variable frozen_value equal ${"Name of original variable that changes in time"}.
##The variable will then be equal to the value of the original compute at the time the new variable is created and the var "frozen_value" will remain constant. (as opposed to v_variable)

##lammps commands diamond,fcc,hcp all construct in the 
##conventional sense (simple cubic lattice)
#lattice diamond ${a_o} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1  
#region box block 0 1 0 1 0 1 units lattice
#create_box   1 box
#create_atoms 1 box
#replicate ${rep_x} ${rep_y} ${rep_z}

###SOME QUANTITIES OF INTEREST
#variable a equal (lx)/(${rep_x}) #INSTANTAEOUS LATTICE CONSTANT 
#compute 1 all pe/atom
#compute 2 all reduce ave c_1

####load potential
#pair_style sw
#pair_coeff * * pot/Si.sw Si1.00194017060324
#mass * 28.0855

###RELAX STRUCTURE
#fix 1 all box/relax iso 0.0 vmax 0.001
#thermo 1 
#thermo_style custom step pe c_2 press pxx pyy pzz v_a
#min_style cg 
#minimize 1e-25 1e-25 5000 10000 

#write_data lammps_data.0K.minimized #data file

#print "minimum-lattice constants: ${a}"

###DEFINE EQ 0K BLOCK SIZE FOR THERMAL EXPANSION CALCULATION
#variable lxo equal ${rep_x}*${a}  #treated as constants 
#variable lyo equal ${rep_y}*${a}
#variable lzo equal ${rep_z}*${a}
#clear

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
#pair_style sw1.00194017060324
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

####--------------------------------------------------------------------
####SIMULATION-3: Scale 0K structure to relevant temperature and
####              randomize positions in NVT (also check pressure state)
####--------------------------------------------------------------------

#ASSUMING SCALE FACTOR FOR CORRECT TEMP IS KNOWN YOU CAN SKIP TO HERE
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
fix ave_press all ave/time 1 ${N_PROD_NVT} ${N_PROD_NVT} v_Pxx v_Pyy v_Pzz v_Pxy v_Pxz v_Pyz  file  ${LOC1}/NVT-ave-pressure-5-5.dat 

run ${N_PROD_NVT} 

write_data ${LOC1}/NVT-EQ-STRUCTURE-5-5 #data file

clear 

####--------------------------------------------------------------------
####SIMULATION-4: RUN NVE AND APPLY MULLER-PLATHE TO GET THERMAL CONDUCTIVITY
####--------------------------------------------------------------------

log ${LOC2}/log-TC-eq-5-5  #save TD output to this file	

units 		metal
dimension	3
boundary	p	p	p
atom_style	atomic

##READ ATOMS
read_data ${LOC1}/NVT-EQ-STRUCTURE-5-5
velocity all create ${T} 4928459 rot yes mom yes dist gaussian

##load potential
pair_style sw
pair_coeff * * pot/Si.sw Si
mass * 28.0855

reset_timestep 0

#LET THE FIX EQUILIBIRATE 
###fix ID group-ID thermal/conductivity N edim Nbin keyword value ...
fix 1 all thermal/conductivity ${swap_every} y ${TC_num_layers} swap ${num_swap}  

fix       2 all nve

thermo ${thermo_save}
thermo_style custom time temp f_1 #v_flux
run ${N_TC_EQ} 

write_data ${LOC2}/TC-EQ-STRUCTURE-5-5 #data file
clear 

####--------------------------------------------------------------------
####SIMULATION-4: RUN NVE AND APPLY MULLER-PLATHE TO GET THERMAL CONDUCTIVITY
####--------------------------------------------------------------------

log ${LOC2}/log-TC-prod-5-5  #save TD output to this file	

units 		metal
dimension	3
boundary	p	p	p
atom_style	atomic

###READ ATOMS
read_data ${LOC2}/TC-EQ-STRUCTURE-5-5

###load potential
pair_style sw
pair_coeff * * pot/Si.sw Si
mass * 28.0855


####fix ID group-ID thermal/conductivity N edim Nbin keyword value ...
reset_timestep 0

fix 1 all thermal/conductivity ${swap_every} y ${TC_num_layers} swap ${num_swap}  

fix       2 all nve

#TEMPERATURE PROFILE
compute   KE all ke/atom
variable  temp1 atom c_KE/(1.5*0.00008617)
compute   layers all chunk/atom bin/1d y lower 2.0 #units reduced
fix 3 all ave/chunk 1 ${N_TC} ${N_TC} layers v_temp1 file ${LOC2}/TC-TEMP-PROFILE-5-5

thermo ${thermo_save}

variable A equal lx*lz
variable flux equal f_1/((time-0.000000000001)*v_A) #ADD -0.000001 to avoid division by zero
thermo_style custom time temp v_A f_1 v_flux

run ${N_TC} 

print ${flux} file ${LOC2}/flux-5-5.dat


####-------------------SOME INFO ON VARIOUS FIXES---------------------------


####write_data snap/data.${N_steps}
#####dump dump1 all custom ${snap_save} snap/snap.*.lammps id type x y z 


##### CONSTRUCT A TEMPERATURE PROFILE 
##### THESE BINS ARE NOT RELATED TO thermal/conductivity Bins
#####  bin/1d args = dim origin delta
#####    dim = x or y or z
#####    origin = lower or center or upper or coordinate value (distance units)
#####    delta = thickness of spatial bins in dim (distance units)



#####fix ID group-ID ave/chunk Nevery Nrepeat Nfreq chunkID value1 value2 ... keyword args ...
#####ID, group-ID are documented in fix command
#####ave/chunk = style name of this fix command
#####Nevery = use input values every this many timesteps
#####Nrepeat = # of times to use input values for calculating averages
#####Nfreq = calculate averages every this many timesteps
#####chunkID = ID of compute chunk/atom command
#####one or more input values can be listed
#####value = vx, vy, vz, fx, fy, fz, density/mass, density/number, temp, c_ID, c_ID[I], f_ID, f_ID[I], v_name

#####For example, if Nevery=2, Nrepeat=6, and Nfreq=100, then values on timesteps 90,92,94,96,98,100 will be used to compute the final aver
####age on timestep 100. Similarly for timesteps 190,192,194,196,198,200 on timestep 200, etc. If Nrepeat=1 and Nfreq = 100, then no time a
####veraging is done; values are simply generated on timesteps 100,200,etc.





#####fix ID group-ID thermal/conductivity N edim Nbin keyword value ...
#####ID, group-ID are documented in fix command
#####thermal/conductivity = style name of this fix command
#####N = perform kinetic energy exchange every N steps
#####edim = x or y or z = direction of kinetic energy transfer
#####Nbin = # of layers in edim direction (must be even number)
#####zero or more keyword/value pairs may be appended
#####keyword = swap
#####swap value = Nswap = number of swaps to perform every N steps


