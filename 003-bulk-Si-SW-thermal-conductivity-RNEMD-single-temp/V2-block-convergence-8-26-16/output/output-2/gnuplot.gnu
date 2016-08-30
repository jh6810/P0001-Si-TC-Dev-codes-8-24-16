f(x)=A1*x+A0
fit f(x)'y-vs-T.dat' u 1:2 via A1,A0
set print "slope.dat" append
print A1 
#K/angstrom

set term png   
set output "plot.png"
plot f(x), 'y-vs-T.dat' u 1:2 pt 7 ps 1.5 , 'temp.dat' u 2:4 


