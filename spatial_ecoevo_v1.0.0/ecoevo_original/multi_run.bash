#!/bin/bash

seed_array=()
run_name="PeriodicCC"
y_array=()
x_array=()
Cmax=()
for ((i=1; i<=4; i++)); do
   y_array[i]=100
   x_array[i]=1
   seed_array[i]=$((3695*i+2**i))
   Cmax[i]=$((4*i+7))
   Cmin[i]=$((2*i+3))
done

for ((i=1; i<=2; i++)); do
  for ((c=2; c<=5; c++)); do
    for ((t=1; t<=4; t++)); do
  	./run.bash ${seed_array[$i]} "$run_name${seed_array[$i]}_AbsSin_${Cmax[$t]}_${Cmin[$t]}" ${y_array[$i]} ${x_array[$i]} $c "FALSE" ${Cmax[$t]} ${Cmin[$t]} 
  	./run.bash ${seed_array[$i]} "$run_name${seed_array[$i]}_FullSin_${Cmax[$t]}_${Cmin[$t]}" ${y_array[$i]} ${x_array[$i]} $c "TRUE" ${Cmax[$t]} ${Cmin[$t]}
	  done
	done
done


