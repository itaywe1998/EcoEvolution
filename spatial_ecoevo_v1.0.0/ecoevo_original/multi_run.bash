#!/bin/bash

seed=1234
run_name="MultiDualRun"
y_array=()
x_array=()
for ((i=1; i<=1; i++)); do
   y_array[i]=100
   x_array[i]=$((5**i))
done

for ((i=1; i<=1; i++)); do
	./run.bash $seed $run_name ${y_array[$i]} ${x_array[$i]}
done


