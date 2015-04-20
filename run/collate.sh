#! /bin/bash

RESFILE=$5
REPS=$(( $2/$3 ))
FOLDERS=$3
RUNF="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"   # run folder
NUMBER=$4
FILE="MCE$NUMBER.sh"
EXDIR=$1
out=${EXDIR##*/}
if [[ ! -d "$EXDIR" ]]; then
 echo "Error. Execution folder does not exist. Check that you are running the correct instance of the result shell script."
 exit 1
fi

folseq=( `seq 1 1 $FOLDERS` )
if [[ ! -z $( command -v qstat ) ]]; then HSTFLG=1; else HSTFLG=0; fi    # HSTFLG=0 ---> qstat is present
cd ..
OUTF1="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESDIR="$OUTF1/$out-$2-reps-$FOLDERS-parts"
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
  if [[ -f Clonetrack-001.out ]]; then
   for j in Clonetrack*; do
    mv $j $RESDIR/${j%.out}_$i.out
   done
  fi
  grep -i "prop NO" input.dat > /dev/null
  propchk=$?
  if [[ propchk -eq 0 ]]; then
   for p in Outbs-*.out ; do mv $p "$RESDIR"/${p%.out}_$i.out; done
  fi
  if [[ -f $SUBDIR/Outbs-001-00000-0.out ]]; then
   echo "yes"
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
 cp $RUNF/timehist.exe .
 ./timehist.exe $FOLDERS $2
fi

if [[ -f normpop_1.out ]]; then
 cols=$( awk '{print NF}' normpop_1.out | sort -nu | tail -n 1 )
 cp $RUNF/avrgpops.exe .
 ./avrgpops.exe $FOLDERS $2 $cols
fi

for i in "${folseq[@]}"; do
 SUBDIR="$EXDIR/$i-run"
 if [[ ! -d "$SUBDIR" ]]; then    #if directory doesn't exist
  echo "Error. Expected folder $SUBDIR does not exist here"
  exit 1
 else
  rm -rf $SUBDIR
 fi
done

GPL=1 
command -v gnuplot >/dev/null 2>&1 || { echo >&2 "Unable to plot results as gnuplot is missing."; GPL=0; }
cd "$RESDIR"
if [[ $GPL -eq 1 ]]; then
 for i in *.gpl; do gnuplot $i; done
fi

cd $RESDIR
if [[ -f "avrgpops.exe" ]]; then rm avrgpops.exe; fi
if [[ -f "timehist.exe" ]]; then rm timehist.exe; fi
if [[ $HSTFLG -eq 1 ]]; then
 cd $EXDIR
 for j in $FILE.[oe]*; do mv $j $RESDIR/$j.out; done
fi
cd "$RUNF"
rm -rf $EXDIR
if [[ -f "interpolate.exe" ]]; then rm interpolate.exe; fi
if [[ -f "MCESB.exe" ]]; then rm MCESB.exe; fi
if [[ -f "$RESFILE" ]]; then rm $RESFILE; fi
if [[ -f "$FILE" ]]; then rm $FILE; fi
if [[ -f "$FILE.[oe]*" ]]; then
 for i in "$FILE.[oe]*" ]]; do mv $i $OUTF1/; done
fi 
