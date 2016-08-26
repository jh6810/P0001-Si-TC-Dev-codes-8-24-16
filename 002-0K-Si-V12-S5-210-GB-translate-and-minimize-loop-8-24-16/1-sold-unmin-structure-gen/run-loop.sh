rm -rf runs 
mkdir runs 
rm -rf structures-plt
mkdir structures-plt
rm -rf structures-lammps
mkdir structures-lammps

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
echo $x\_$z
cd runs
rm -rf $x\_$z
cp -r ../template $x\_$z
cd $x\_$z
	sed "s/x_in/$x/g" defects.def  >tmp
	sed "s/z_in/$z/g" tmp  >tmp2
	rm tmp
	mv tmp2 defects.def
	more defects.def
	./run.sh
	cp out/config.plt ../../structures-plt/$x\_$z.plt
	mv out/config.plt ../../plt2lammps-07-18-13/temp.plt
cd ../
cd ../
cd plt2lammps-07-18-13
./a.out -i temp.plt -t Si > temp.lammps
mv temp.lammps ../structures-lammps/$x\_$z.lammps
rm temp.plt
cd ../

done
done 
