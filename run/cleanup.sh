#! /bin/bash

REPS=$(( $1/$2 ))
FOLDERS=$2
RUNF="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"   # run folder
NUMBER=$3
FILE="MCE$NUMBER.sh"
cd ..
#if [[ -d Reset ]]; then
# cd Reset
# rm -rf *
# cd ..
#else
# mkdir Reset
#fi
#cp -rf EXEC/ Reset/
cd run
#cp avrgpops.exe exdir.dat interpolate.exe $FILE MCESB.exe result.sh runconds.dat timehist.exe ../Reset
folseq=( `seq 1 1 $FOLDERS` )
if [[ ! -z $( command -v qstat ) ]]; then HSTFLG=1; else HSTFLG=0; fi
cd ..
OUTF1="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTF="$OUTF1/MCE-$1-reps-$FOLDERS-parts-$NUMBER"
if [[ ! -d $OUTF ]]; then mkdir $OUTF; fi
cd $RUNF
if [[ -f exdir.dat ]]; then EXDIR1=$( head -n1 exdir.dat ); else echo "exdir.dat is missing!"; exit 1; fi
for a in 1 2; do
 cd $RUNF
 RESDIR="$OUTF/output-MCEv$a"
 EXDIR="$EXDIR1/MCEv$a-$NUMBER"
 if [[ ! -d "$EXDIR" ]]; then
  echo "Error. Execution folder does not exist. Check that number is correct."
  exit 1
 fi
 if [[ ! -d "$RESDIR" ]]; then    #if directory doesn't exist
  mkdir "$RESDIR"
 else
  cd "$RESDIR"
  read -p "Output Directory Exists! Overwrite old data? (Y=Overwrite, N=Exit) " -n 1 -r
  echo    # move to a new line
  if [[ $REPLY =~ ^[Yy]$ ]]; then
   if [[ "$(ls -A )" ]]; then rm *; fi
  else
   exit 1
  fi
  cd ..
 fi
 for i in "${folseq[@]}"; do
  SUBDIR="$EXDIR/$i-run"
  if [[ ! -d "$SUBDIR" ]]; then    #if directory doesn't exist
   echo "Error. Expected folder $SUBDIR does not exist"
   exit 1
  else
   cd "$SUBDIR"
   if [[ $HSTFLG -eq 0 ]]; then
    if [[ -f $FILE.[oe]* ]]; then
     for j in $FILE.o*; do
      mv $j $RESDIR/$j.out
     done
    fi
   fi
   if [[ -f normpop.out ]]; then
    mv normpop.out "$RESDIR"/normpop_$i.out
   fi
   if [[ -f Outbs-001.out ]]; then
    for p in Outbs-*.out ; do mv $p "$RESDIR"/${p%.out}_$i.out; done
   fi
   if [[ -f timehist.out ]]; then
    mv timehist.out "$RESDIR"/timehist_$i.out
    mv timesteps.out "$RESDIR"/timesteps_$i.out
   fi
  fi
  cd ..
 done
 cd "$RESDIR"
 cp $SUBDIR/input.dat $SUBDIR/inham.dat $SUBDIR/prop.dat .
 if [[ -f timehist_1.out ]]; then 
  for i in "${folseq[@]}"; do cat "timesteps_$i.out"; done > timesteps.out
  for i in "${folseq[@]}"; do rm "timesteps_$i.out"; done
  cp ../../run/timehist.exe .
  ./timehist.exe $FOLDERS $1
 fi
 if [[ -f normpop_1.out ]]; then
  cp ../../run/avrgpops.exe .
  ./avrgpops.exe $FOLDERS $1
 fi
 for i in "${folseq[@]}"; do
  SUBDIR="$EXDIR/$i-run"
  if [[ ! -d "$SUBDIR" ]]; then    #if directory doesn't exist
   echo "Error. Expected folder $SUBDIR does not exist here"
   exit 1
  else
   cd "$SUBDIR"
   if [[ "$(ls -A )" ]]; then rm *; fi
   cd ..
   rmdir "$i-run"
  fi
 done
 GPL=1 
 command -v gnuplot >/dev/null 2>&1 || { echo >&2 "Unable to plot results as gnuplot is missing."; GPL=0; }
 cd "$RESDIR"
 if [[ $GPL -eq 1 ]]; then
  if [[ -f plothist.gpl ]]; then gnuplot plothist.gpl; fi
  if [[ -f plothistall.gpl ]]; then gnuplot plothistall.gpl; fi
  if [[ -f plotpopdiff.gpl ]]; then gnuplot plotpopdiff.gpl; fi
  if [[ -f plotpopres.gpl ]]; then gnuplot plotpopres.gpl; fi
 fi
 cd $RESDIR
 if [[ -f "avrgpops.exe" ]]; then rm avrgpops.exe; fi
 if [[ -f "timehist.exe" ]]; then rm timehist.exe; fi
 if [[ $HSTFLG -eq 1 ]]; then
  cd $EXDIR
  for j in $FILE.[oe]*; do mv $j $RESDIR/$j.out; done
 fi
 cd $EXDIR
 if [[ "$(ls -A )" ]]; then rm *; fi
 cd ..
 rmdir $EXDIR
 cd "$RUNF"
done
cd $RUNF
#if [[ -f "avrgpops.exe" ]]; then rm avrgpops.exe; fi
#if [[ -f "timehist.exe" ]]; then rm timehist.exe; fi
if [[ -f "interpolate.exe" ]]; then rm interpolate.exe; fi
if [[ -f "MCESB.exe" ]]; then rm MCESB.exe; fi
#if [[ -f "exdir.dat" ]]; then rm exdir.dat; fi
if [[ -f "result.sh" ]]; then rm result.sh; fi
if [[ -f "$FILE" ]]; then rm $FILE; fi
if [[ -f "$FILE.[oe]*" ]]; then
 for i in "$FILE.[oe]*" ]]; do mv $i $OUTF1/; done
fi 
#cd $EXDIR1
#if [[ "$(ls -A )" ]]; then rm *; fi
#cd ..
#if [[ -d "EXEC" ]]; then rmdir EXEC; fi 
#cd "$RUNF"
