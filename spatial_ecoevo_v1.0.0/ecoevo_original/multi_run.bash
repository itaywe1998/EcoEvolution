#!/bin/bash

seed=1234
run_name="MultiDualRun"
y_array=()
for ((i=1; i<=5; i++)); do
   y_array[i]=$((i*50))
done

for y in "${y_array[@]}"
do
	./run.bash $seed $run_name $y
done

