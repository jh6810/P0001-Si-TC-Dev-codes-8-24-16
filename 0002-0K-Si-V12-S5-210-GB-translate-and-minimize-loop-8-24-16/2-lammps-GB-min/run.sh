rm -rf logs 
mkdir logs 

rm -rf min-structures 
mkdir min-structures 

rm Tx-Tz-gbejm2.dat
#awk -v min=1000 '{if($3<min){min=$3}{Tx=$1}{Tz=$2}} END {print Tx,Tz,min}' 'Tx-Tz-gbejm2.dat'

#ACTUAL S5 PERIODIC TRANSLATIONS 
#for x in 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0 10.5 11.0 11.5 12.0 12.5
#do 
#for z in 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 
#do

#TEST CASE
for x in 0.0 
do 
for z in 0.0 
do 
	Tx=$x
	Tz=$z
	i=$Tx\_$Tz.lammps 
	cd template 
	cp ../structures-lammps/$i ./structure.lammps
	./run.sh

	En=$(awk '{print $0}' 'gbejm2.dat')
	echo $Tx $Tz $En >> ../Tx-Tz-gbejm2.dat

	cp snap/lammps_data.min ../min-structures/$Tx\_$Tz.lammps
	cp log.lammps ../logs/$Tx\_$Tz.log
	cd ../
done 
done
