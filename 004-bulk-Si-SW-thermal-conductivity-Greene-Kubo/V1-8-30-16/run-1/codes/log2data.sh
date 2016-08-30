
#stop=1000000 #incase simulation is still running
start=$(awk '{if($1 =="Time"){{line=NR} {print line}}}' 'snap/L1.log')
stop=$(awk '{if($1 =="Loop"){{line=NR} {print line}}}' 'snap/L1.log')

awk -v start=$start -v stop=$stop '{if(NR>start && NR<stop ){print $0}}' 'snap/L1.log' > snap/TD.dat 

