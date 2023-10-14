#!/bin/bash

seed_array=()
run_name="BiggerCompetitionDiffSeeds"
y_array=()
x_array=()
for ((i=1; i<=7; i++)); do
   y_array[i]=100
   x_array[i]=1
   seed_array[i]=$((721*i+3**i))
done

for ((i=1; i<=7; i++)); do
	./run.bash ${seed_array[$i]} $run_name ${y_array[$i]} ${x_array[$i]}
done


