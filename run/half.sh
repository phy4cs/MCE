#! /bin/bash

#######There are now no parameters in this file which need manually changing. Everything that needs changing between runs is in the input.dat file

REPS=$(( $1/$3 ))     # number of repeats per folder
CORES=$2
FOLDERS=$3
RUNF=$PWD   # run folder
if [[ -d ../build ]]; then
 cd ../build
 BUILD=$PWD
 cd $RUNF
else
 echo "No build folder! What's going on? I'm scared."
 echo ""
 echo ""
 echo ""
 echo "Hold me."
 exit 1
fi
NUMBER=$RANDOM
FILE="MCE$NUMBER.sh"
if [[ ! -z $( command -v qstat ) ]]; then HSTFLG=1; else HSTFLG=0; fi
if [[ -n $( echo $HOSTNAME | fgrep -e "arc1" -e "polaris" -e "arc2" ) ]]; then HPCFLG=1; else HPCFLG=0; fi
if [[ $HPCFLG -eq 0 ]]; then
 cd ..
 if [[ ! -d "EXEC" ]]; then
  mkdir EXEC
  cd EXEC
 else
  cd EXEC
 fi 
 EXDIR1=$PWD    # execution folder
else
 EXDIR1="/nobackup/$LOGNAME"
fi
if [[ ! -d "$EXDIR1" ]]; then echo "Cannot find execution directory $EXDIR1. Exitting"; exit 1; fi
folseq=( `seq 1 $FOLDERS` )
cd $BUILD
if [[ HPCFLG -eq 1 ]]; then
 make -f makefile_arc
else
 make -f makefile_chmlin
fi 
cp *.exe $RUNF 
cd $RUNF
if [[ $? -ne 0 ]]; then 
 echo "Compilation Error! Exitting"
 exit 1
fi
method=`grep -i "^method" input.dat`
if [[ $? == 1 || $? == 2 ]]; then
 echo "Could not read the method from input.dat. Exitting"
 exit 1
fi
method=${method#* }
if [[ $method == "MCE12" ]]; then
   k=2
else
   k=1
fi
sed -i "s/^Repeats.*/Repeats $REPS/g" input.dat
grep "^Repeats $REPS" input.dat > /dev/null
if [[ $? == 1 || $? == 2 ]]; then
 echo "Could not change the number of repeats in input.dat. Exitting"
 exit 1
fi 
methseq=( `seq 1 $k` )
outfol=`grep -i "^Runfolder" input.dat`
if [[ $? == 1 || $? == 2 ]]; then
 echo "Could not read the execution folder from input.dat. Exitting"
 exit 1
fi
outfol=${outfol#* }
sys=`grep -i "^System:" input.dat`
if [[ $? == 1 || $? == 2 ]]; then
 echo "Could not read the system from input.dat. Exitting"
 exit 1
fi
sys=${sys#* }
echo $outfol | grep -i "Default" > /dev/null
outdef=$?
for a in "${methseq[@]}"; do
 if [[ $k == 2 ]]; then
  sed -i "s/^method.*/method MCEv$a/g" input.dat
  if [[ $outdef == 0 ]]; then
   outfol2="MCEv$a-$sys-$NUMBER"
  else
   outfol2="MCEv$a-$sys-$outfol"
  fi
 else
  if [[ $outdef == 0 ]]; then
   outfol2="$method-$sys-$NUMBER"
  else
   outfol2="$method-$sys-$outfol"
  fi
 fi
 EXDIR="$EXDIR1/$outfol2"
 if [[ ! -d $EXDIR ]]; then mkdir $EXDIR; fi
 for i in "${folseq[@]}"; do
  SUBDIR="$EXDIR/$i-run"
  if [[ ! -d "$SUBDIR" ]]; then    #if directory doesn't exist
   mkdir "$SUBDIR"
  else
   cd "$SUBDIR"
   if [[ "$(ls -A )" ]]; then rm *.*; fi             #remove old run files
  fi
 done
 cd "$RUNF"
 chk=`grep "gen YES" input.dat`
 if [[ ! -z ${chk} ]]; then 
  gen=1
 else
  gen=0
 fi
 cd ../run/
 echo "#$ -cwd -V -m e" > $FILE
 if [[ $CORES -ne 1 ]]; then echo "#$ -pe smp $CORES" >> $FILE; fi
 echo "#$ -l h_rt=40:00:00" >> $FILE
 echo "#$ -l h_vmem=4G" >> $FILE
 echo "#$ -t 1-$FOLDERS" >> $FILE
 echo "cd $EXDIR/"'$SGE_TASK_ID'"-run/" >> $FILE
 echo "echo "'"Running on $HOSTNAME in folder $PWD"' >> $FILE 
 if [[ $HPCFLG -eq 1 ]]; then
  echo "module load mkl" >> $FILE
 fi
 echo "./MCESB.exe" >> $FILE
 for i in "${folseq[@]}"; do
  SUBDIR="$EXDIR/$i-run"
  cd "$RUNF"
  cp inham.dat input.dat MCESB.exe prop.dat $SUBDIR/
  if [[ $gen -eq 0 ]]; then
   if [[ -f "Outbs-001_$i.out" ]]; then 
    echo "Outbs-001_$i.out found in $PWD"
    for x in Outbs-*_$i.out; do
     cp $x $SUBDIR/${x%_$i.out}.out
    done
   else
    echo "Outbs-001_$i.out not found in $PWD"
    echo "For propagation to occur without basis set generation, all relevant input bases must be present"
    exit 1
   fi
  fi
  if [ $HSTFLG -eq 0 ]; then
   cd $SUBDIR/
   echo "Program Executing in $EXDIR"
   ./MCESB.exe #&> $FILE.o1 &
   cd $RUNF
  fi
 done
 mv $FILE $EXDIR
 cd $EXDIR
 if [[ $CORES -ne 1 ]]; then export OMP_NUM_THREADS=$CORES; fi
 if [[ $HSTFLG -eq 1 ]]; then
  qsub $FILE
 fi
 cd "$RUNF"
done
echo "./cleanup.sh $1 $3 $NUMBER" > result.sh
if [[ $k == 2 ]]; then
 sed -i "s/^method.*/method MCE12/g" input.dat
fi
chmod u+x result.sh
