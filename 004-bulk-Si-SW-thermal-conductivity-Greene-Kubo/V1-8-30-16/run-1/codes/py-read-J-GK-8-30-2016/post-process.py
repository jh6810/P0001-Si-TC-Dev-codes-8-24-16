import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq,rfft,rfftfreq,ifft
import random
import os 
import shutil
from pylab import *
import pylab
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import time
#------------------------------GET VOLUME------------------------------
file_loc='../../output/output-2/volume.dat'

input_file = open( file_loc, "r" )

#STRIP AWAY BLANK LINES
with open( file_loc, "r" ) as f:
	for line in f: 
		Vol=float(line)
print Vol,"angstrom^3"
#------------------------------READ LOG FILE ------------------------------
file_loc='../../output/output-2/NVE.log'

input_file = open( file_loc, "r" )

#STRIP AWAY BLANK LINES
with open( file_loc, "r" ) as f:
    names_list = [line.strip() for line in f if line.strip()]

#LAMMPS VARIABLE NAMES AND START AND STOP LINE
line_num=1
for i in range(0,len(names_list)):
	parts = names_list[i].split()
	if(parts[0]=='Time'): 
		line_1=line_num
		print "SET OF ALL LAMMPS VARIABLES:"
		print parts


		keyline=[]
		#READ SELECTED LAMMPS DATA
		list_1=[0,1,2,8] 
		for j in range(0,len(list_1)):
			keyline.append(str(parts[list_1[j]]))
			exec(str(parts[list_1[j]]) + " = []")


	if(parts[0]=='Loop'): 
		line_2=line_num

#	IF JOB IS STILL RUNNING
#	line_2=line_num-1

	line_num=line_num+1

#READ LESS OF THE FILE IF DESIRED 
#line_2=int(line_2/10)

#print "BEGIN READING LAMMPS VARIABLES:"
print keyline

###READ ALL DATA INTO RESPECTIVE ARRAYS 
load_command=str(keyline[0]) + ".append(float(parts["+str(list_1[0])+"]));"
for j in range(1,len(keyline)):
	load_command=load_command+str(keyline[j]) + ".append(float(parts["+str(list_1[j])+"]));"
i=line_1+1
while(i<line_2-1): 
	parts = names_list[i].split()
	exec(load_command)
	i=i+1;
print "DONE READING:"

##---------------------------ENSEMBLE AVERAGES---------------------------

####SAVE OUTPUT
if os.path.exists('output'):
    shutil.rmtree('output')
os.makedirs('output')
os.makedirs('output/plots')

t=np.array(Time); 

#MEAN IS ALREADY BASICALLY ZERO (BUT SUBTRACT OUT ANYWAY)

J=np.array(J)#-np.mean(J); 
KE=np.array(KinEng)-np.mean(KinEng); 

T=np.mean(Temp); print T

##BREAK DATA SET INTO CHUNKS TO TAKE ENSEMBLE AVERAGE
dt1=t[2]-t[1]
to=t[0]
tf=t[len(t)-1]
Np=len(t)
dt_window=100 #in ps  
N_windows=int((tf-to)/dt_window)

di=int(floor(Np/N_windows)) #needs to be odd
if(di%2==0):
	di=di-1

print "Number of windows:",N_windows
print "length of window (ps):",dt_window
print "points per window:",di

t_supper=[]; j_ACF_supper=[]; 
f_supper=[]; j_SPD_supper=[];

k=0; N2=0.0
while(k<Np): 
	tsub=np.array(t[k:k+di])
	KEsub=np.array(KE[k:k+di])

	Jsub=np.array(J[k:k+di])

	if(len(tsub)==di and len(Jsub)==di):
		print "len(tsub),len(Jsub)",len(tsub),len(Jsub)
		tsub=np.array(tsub)-tsub[0]
		print "tsub_f",tsub[len(tsub)-1]


		N=len(tsub); dt=tsub[2]-tsub[1]; f=fftfreq(N)/(dt);
		print "dt (ps):",dt

		J_SPD_i=(np.absolute(fft(Jsub)/N)**2.0)
		KE_SPD_i=(np.absolute(fft(KEsub)/N)**2.0)

		J_ACF_i=np.real(N*ifft(J_SPD_i))
		KE_ACF_i=np.real(N*ifft(KE_SPD_i))



		for j in range(0,len(f)):
			f_supper.append(f[j])
			j_SPD_supper.append(J_SPD_i[j])

		for j in range(0,len(tsub)):
			t_supper.append(tsub[j])
			j_ACF_supper.append(J_ACF_i[j])

	#	#-------SAVE KE ACF DATA---------
		if(N2==0):
			t_ave=np.array(tsub)
			f_ave=np.array(f)

			J_SPD_ave=np.array(J_SPD_i)
			KE_SPD_ave=np.array(KE_SPD_i)

			J_ACF_ave=np.array(J_ACF_i)
			KE_ACF_ave=np.array(KE_ACF_i)
		else: 
			f_ave=np.add(f_ave,np.array(f))
			t_ave=np.add(t_ave,np.array(tsub))

			J_SPD_ave=np.add(J_SPD_ave,np.array(J_SPD_i))
			KE_SPD_ave=np.add(KE_SPD_ave,np.array(KE_SPD_i))

			J_ACF_ave=np.add(J_ACF_ave,np.array(J_ACF_i))
			KE_ACF_ave=np.add(KE_ACF_ave,np.array(KE_ACF_i))

		N2=N2+1.0

	k=k+di

t_ave=np.array(t_ave/N2); 
J_ACF_ave=np.real(J_ACF_ave/N2); KE_ACF_ave=np.array(KE_ACF_ave/N2);
f_ave=np.array(f_ave/N2); 
J_SPD_ave=np.array(J_SPD_ave/N2); KE_SPD_ave=np.array(KE_SPD_ave/N2);

##SOME PLOTS
#plt.plot(f_supper,j_SPD_supper,'o',f_ave,J_SPD_ave,'.'); show();  plt.clf(); 
#plt.plot(t_supper,j_ACF_supper,'o',t_ave,J_ACF_ave,'.'); show();  plt.clf(); 

#plt.plot(f_ave,KE_SPD_ave,'.'); show();  plt.clf(); 
#plt.plot(f_ave/2.0,KE_SPD_ave,'.'); show();  plt.clf(); 
#plt.plot(t_ave,KE_ACF_ave,'.'); show();  plt.clf(); 

#1 ev/(ang*ps)=1602.1766 W/m
print "thermal conductivity (W/(m*K)) (expected at 500K):",80
kappa=1602.17662*sum(J_ACF_ave[0:len(J_ACF_ave)-2])*(t_ave[2]-t_ave[1])*Vol/((T**2.0)*8.61733*10**-5)/2.0 #divid by two (periodicity of ACF) 
print "thermal conductivity (W/(m*K))(LHS):",kappa


with open("output/ACF-AVE.dat", "a") as myfile:
	for i in range(0,int(len(t_ave)/2)):
		myfile.write(' %13.6e %13.6e %13.6e \n' % (t_ave[i],np.real(J_ACF_ave[i]),np.real(KE_ACF_ave[i])))
myfile.close()

with open("output/SPD-AVE.dat", "a") as myfile:
	for i in range(0,int(len(f_ave))):
		myfile.write(' %13.6e %13.6e %13.6e \n' % (f_ave[i],J_SPD_ave[i],KE_SPD_ave[i]))
myfile.close()


##------------------------------EXTRA CODES------------------------------



##plt.plot(t,(Jx-np.mean(Jx))/np.std(Jx),'-',t,(KinEng-np.mean(KinEng))/np.std(KinEng),'.-',t,(TotEng-np.mean(TotEng))/np.std(TotEng),'.'); show();  plt.clf();
##plt.plot(t,(Jx-np.mean(Jx))/np.std(Jx),'-',t,(KinEng-np.mean(KinEng))/np.std(KinEng),'-'); show();  plt.clf();
##plt.plot(t,(Jx-np.mean(Jx)),'-',t,(KinEng-np.mean(KinEng)),'o',t,(TotEng-np.mean(TotEng)),'.'); show();  plt.clf(); 




###FIGURE OUT #plt.plot(t,(Jx-np.mean(Jx))/np.std(Jx),'-',t,(KinEng-np.mean(KinEng))/np.std(KinEng),'.-',t,(TotEng-np.mean(TotEng))/np.std(TotEng),'.'); show();  plt.clf();
#plt.plot(t,(Jx-np.mean(Jx))/np.std(Jx),'-',t,(KinEng-np.mean(KinEng))/np.std(KinEng),'-'); show();  plt.clf();
#plt.plot(t,(Jx-np.mean(Jx)),'-',t,(KinEng-np.mean(KinEng)),'o',t,(TotEng-np.mean(TotEng)),'.'); show();  plt.clf(); 



###TRAPIZOID RULE
##I=0.0; dt=(t_ave[2]-t_ave[1])
##for i in range(0,len(ACF_ave)-1): #EXCLUSIVE ON UPPER LIMIT
##	I=I+dt*(ACF_ave[i+1]+ACF_ave[i])/2.0



#plt.plot(t,(Jx-np.mean(Jx))/np.std(Jx),'.',t,(J-np.mean(J))/np.std(J),'-'); plt.savefig('output/plots/plot-1.png');  plt.clf();








