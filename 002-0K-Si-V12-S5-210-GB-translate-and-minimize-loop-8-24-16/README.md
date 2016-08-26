# 0002-0K-Si-V12-S5-210-GB-translate-and-minimize-loop-8-24-16
# contact: James Hickman (GMU) (jhickma3@gmu.edu)
# Requirement: you need to have lammps installed on computer (see lammps.com) 
# Requirement: you need to have SOLD installed on computer (contact email above) 

# This code attempts to find the mininmum energy Sigma 5 210 grain boundary for 
# the Si potential included in the code set.This is done by generating many structures
# with a fixed cross sectional area corresponding to the stress free lattice but with 
# translations of one grain relative to another in the in-plane GB directions. These 
# structures are then feed into lammps to be minimized with respect to local atomic
# displacements followed by block deformation in the GB normal direction 
# (minimization method was conjugate gradient)  (structures are generated using SOLD)
