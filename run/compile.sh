#!/bin/bash
set -e
cd ../src/
if [[ -n $( echo $HOSTNAME | fgrep -e "arc1" -e "polaris" -e "arc2" ) ]]; then ARCFLG=1; else ARCFLG=0; fi
if [[ $ARCFLG -eq 1 ]]; then
 echo "Compiling with Intel Compiler"
 FC=ifort
 if [[ $1 -eq 1 ]]; then
  FFLAGS="-O0 -g"
 else
  echo "Compiling with OpenMP functionality"
  FFLAGS="-openmp -O2"
 fi
else
 FC=gfortran
 echo "Compiling with GNU Compiler"
 if [[ $1 -eq 1 ]]; then
  FFLAGS="-O0 -g"
 else
  echo "Compiling with OpenMP functionality"
  FFLAGS="-fopenmp -O2"
 fi
fi

if [[ $ARCFLG -eq 0 ]]; then
 $FC $FFLAGS -c BLAS.f
 $FC $FFLAGS -c zgesv_with_dpndnss.f
 $FC $FFLAGS -c zheev_with_dpndnss.f
fi
$FC $FFLAGS -c vars.f90
$FC $FFLAGS -c dcerf.f90
$FC $FFLAGS -c sb.f90
$FC $FFLAGS -c hp.f90
$FC $FFLAGS -c fp.f90
$FC $FFLAGS -c mp.f90
$FC $FFLAGS -c iv.f90
$FC $FFLAGS -c cp.f90
$FC $FFLAGS -c hh.f90
$FC $FFLAGS -c redirect.f90
$FC $FFLAGS -c alarrays.f90
$FC $FFLAGS -c outputs.f90
$FC $FFLAGS -c Ham.f90
$FC $FFLAGS -c Chks.f90
$FC $FFLAGS -c randgen.f
$FC $FFLAGS -c readpars.f90
$FC $FFLAGS -c derivsMCE.f90
$FC $FFLAGS -c propMCE.f90
$FC $FFLAGS -c bsetgen.f90
$FC $FFLAGS -c bsetalter.f90
$FC $FFLAGS -c MainMCE.f90

if [[ $ARCFLG -eq 1 ]]; then
 $FC $FFLAGS -o MCESB.exe vars.o dcerf.o sb.o hp.o fp.o mp.o iv.o cp.o hh.o redirect.o alarrays.o outputs.o Ham.o Chks.o randgen.o readpars.o derivsMCE.o propMCE.o bsetgen.o bsetalter.o MainMCE.o -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
else
 $FC $FFLAGS -o MCESB.exe BLAS.o zgesv_with_dpndnss.o zheev_with_dpndnss.o vars.o dcerf.o sb.o hp.o fp.o mp.o iv.o cp.o hh.o redirect.o alarrays.o outputs.o Ham.o Chks.o randgen.o  readpars.o derivsMCE.o propMCE.o bsetgen.o bsetalter.o MainMCE.o -lm -Warray-bounds
fi

$FC $FFLAGS -o avrgpops.exe avrgpops.f90
$FC $FFLAGS -o timehist.exe timehist.f90
#$FC $FFLAGS -o interpolate.exe neville.f90 interpolate.f90

rm *.mod *.o

mv MCESB.exe ../run/
mv avrgpops.exe ../run/
mv timehist.exe ../run/
#mv interpolate.exe ../run/
cd ../run/
