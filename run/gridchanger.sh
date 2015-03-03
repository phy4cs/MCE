#! /bin/bash

qsize=140
psize=35
i=0

if [ ! -d calibinputs ]; then mkdir calibinputs; fi

for grsp in `seq 1.20 0.01 1.50`; do
 dim1=$(perl -w -e "use POSIX; print ceil($qsize/$grsp), qq{\n}")
 dim2=$(perl -w -e "use POSIX; print ceil($psize/$grsp), qq{\n}")
 if [ $(( $dim1 % 2 )) == 1 ]; then dim1=$[$dim1+1]; fi
 if [ $(( $dim2 % 2 )) == 1 ]; then dim2=$[$dim2+1]; fi
 dim3=$(( $dim1 * $dim2 ))
 sed -i "s/^qsizez.*/qsizez $dim1/g" input2.dat
 sed -i "s/^psizez.*/psizez $dim2/g" input2.dat
 sed -i "s/^gridsp.*/gridsp $grsp/g" input2.dat
 sed -i "s/^in_nbf.*/in_nbf $dim3/g" input2.dat 
 for j in `seq 25 25 1000`; do
  sed -i "s/^ALCMP.*/ALCMP $j/g" input2.dat
  sed -i "s/^Runfolder.*/Runfolder ${dim1}x${dim2}-$grsp-$j/g" input2.dat
   i=$[$i+1]
  cp input2.dat ./calibinputs/input.$i
 done
done 
