#!/bin/bash

seed_array=()
run_name="PeriodicCC_Factory"
y_array=()
x_array=()
Cmax=()
for ((i=1; i<=8; i++)); do
   y_array[i]=100
   x_array[i]=1
   seed_array[i]=$((3696*i+2**i))
   Cmax[i]=$((30))
   Cmin[i]=$((30))
done

for ((i=2; i<=2; i++)); do
  for ((c=8; c<=8; c++)); do
  	./run.bash ${seed_array[$i]} "$run_name"_cycles"$c"_seed"${seed_array[$i]}"_AbsSin_"${Cmax[$i]}"_"${Cmin[$i]}" ${y_array[$i]} ${x_array[$i]} $c "FALSE" ${Cmax[$i]} ${Cmin[$i]} 
  #	./run.bash ${seed_array[$i]} "$run_name"_cycles"$c"_seed"${seed_array[$i]}"_FullSin_"${Cmax[$i]}"_"${Cmin[$i]}" ${y_array[$i]} ${x_array[$i]} $c "TRUE" ${Cmax[$i]} ${Cmin[$i]}
	done
done


