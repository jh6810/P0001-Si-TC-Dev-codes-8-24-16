f(x)=A1*x+A0
fit f(x)'y-vs-T.dat' u 1:2 via A1,A0
set print "slope.dat" append
plot f(x), 'y-vs-T.dat' u 1:2 
pause 5
#K/angstrom
print A1 
