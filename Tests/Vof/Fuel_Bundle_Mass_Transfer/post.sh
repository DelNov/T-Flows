gnuplot volume.gnu
gnuplot massTransfer.gnu

cat out |grep Courant > courant.out
sed 's/                       #    //' courant.out > tmp.out
mv -f tmp.out courant.out
gnuplot courant.gnu

