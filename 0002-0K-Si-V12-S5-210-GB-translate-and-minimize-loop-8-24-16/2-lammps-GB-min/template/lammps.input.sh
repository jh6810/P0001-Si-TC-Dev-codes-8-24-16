units 		metal
dimension	3
boundary	p	p	p
atom_style	atomic


variable N_steps equal 1000 
variable thermo_out equal 1 

variable minimumenergy equal -4.630371

compute eng all pe/atom 
compute eatoms all reduce sum c_eng 

read_data structure.lammps

##load potential
pair_style tersoff/mod
pair_coeff * * pot/Si.tersoff-V12.mod Si 
mass * 28.0855


# ---------- Run Minimization --------------------- 
reset_timestep 0 
thermo 10 
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms 
min_style cg 
minimize 1e-15 1e-15 5000 5000 


dump dump1 all custom 50 snap/snap.*.lammps id type x y z 

# ---------- Run Minimization 2--------------------- 
# Now allow the box to expand/contract perpendicular to the grain boundary
reset_timestep 0 
thermo 10 
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms 
fix 1 all box/relax y 0 vmax 0.001
min_style cg 
minimize 1e-15 1e-15 5000 5000 




# ---------- Calculate GB Energy --------------------- 

variable esum equal "v_minimumenergy * count(all)" 
variable xseng equal "c_eatoms - (v_minimumenergy * count(all))" 
variable gbarea equal "lx * lz * 2"  #2 is for the 2 GB 
variable gbe equal "(c_eatoms - (v_minimumenergy * count(all)))/v_gbarea" 
variable gbemJm2 equal ${gbe}*16021.7733 
variable gbernd equal round(${gbemJm2})
variable gbeJm2 equal ${gbe}*16021.7733/1000.0  
print "GB energy is ${gbemJm2} mJ/m^2" 
print "GB energy is ${gbeJm2} mJ/m^2" 

print ${gbeJm2} file gbejm2.dat





write_data snap/lammps_data.min #data file


