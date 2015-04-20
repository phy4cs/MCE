#! /bin/bash

if [[ ! -f folderlist.dat ]]; then
  echo "Folder list is missing! Aborting"
  exit 1
fi
runf=$PWD

IFS=$'\n' read -d '' -r -a lines < folderlist.dat
nelements=${#lines[@]}
max_index=$[$nelements-1]

outf1=${lines[0]}
outf2=${outf1##*/}
cd ../
if [ -d $outf2 ]; then
  echo "The temporary folder for concatenating the normpop files $outf2 already exists! Aborting"
  exit 1
else
  mkdir $outf2
  cd $outf2
  outf=$PWD
  cd ../run
fi

j=0
for i in ${lines[@]}; do
  cd $i
  folnum[j]=$( ls -d */ | wc -l )
  (( j++ ))
done

j=0
for i in ${folnum[@]}; do
  if [ $i -ne ${folnum[0]} ]; then
    echo "Error! Expected to find ${folnum[0]} subdirectories in ${lines[j]}, but got $i"
    exit 1
  fi
  (( j++ ))
done

folders=${folnum[0]}

cd ${lines[0]}/1-run
num1=$( grep "^Repeats" input.dat )
num2=${num1#* }
steptmp=$( grep "^step" prop.dat )
step=${steptmp#* }
endtmp=$( grep "^time_end" prop.dat )
end=${endtmp#* }
end2=$( echo $end | sed -e 's/d/e/g' )
end=$( echo $end2 | awk '{ printf "%#9.8E", $1 }' )

for i in `seq $folders`; do
  for k in `seq -f "%03g" $num2`; do
    p=1
    for j in ${lines[@]}; do
      cd $j/$i-run
      if [ -f normpop-$k.out ]; then
        cp normpop-$k.out $outf/normpop-${k}_${i}_${p}.out
        (( p++ ))
      fi
    done
  done
done

cd $outf
for k in `seq -f "%03g" $num2`; do
  for i in `seq $folders`; do
    for p in normpop-${k}_${i}_*.out; do
      q1=${p%.out}
      q=${q1##*_}
      if [ $q -eq 1 ]; then
        cp $p normpop-${k}_${i}.out
      else
        tail -n +5 -q $p >> normpop-${k}_${i}.out
      fi
    done
  done
done

for k in `seq -f "%03g" $num2`; do
  for i in `seq $folders`; do
    tmp=$( tail -n 1 normpop-${k}_${i}.out )
    tmp2=${tmp#*  }
    tmp=${tmp2%% *}
    endtmp=$( echo $tmp | awk '{ printf "%#9.8E", $1 }' )
    if [[ $endtmp != $end ]]; then
      echo "The file normpop-${k}_${i}.out does not finish at the expected end time."
      echo "This file will be deleted"
      rm normpop-${k}_${i}.out
    fi
  done
done


n=0
for k in `seq -f "%03g" $num2`; do
  for i in `seq $folders`; do
    if [ -f normpop-${k}_${i}.out ]; then
      n=$[$n+1]
      if [ $n -eq 1 ]; then
        cols=$( awk '{print NF}' normpop-${k}_${i}.out | sort -nu | tail -n 1 )
      fi
      mv normpop-${k}_${i}.out ${lines[max_index]}/$i-run/normpop-${k}.out
    fi
  done
done

cd $runf
rm -rf $outf

path=${lines[max_index]}
NUMBER=${path##*-}
reps=$[$num2*$folders]

for i in `seq $folders`; do
  if [ $step == "adaptive" ]; then
    cp interpolate.exe $path/$i-run
    cd $path/$i-run
    ./interpolate.exe $num2 $cols
    cd $runf
  else
    cp subavrg.exe $path/$i-run
    cd $path/$i-run
    ./subavrg.exe $num2 $cols
    cd $runf
  fi
done
  
echo "./collate.sh $path $reps $folders $NUMBER "'$0' > result.sh
chmod u+x result.sh
./result.sh
