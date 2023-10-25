#!/bin/bash

seed_array=()
run_name="PeriodicCC_Factory"
vbar=$((3/(10**5))) #sort of scientific writing (not supported directly in bash)
dbar=$((1/(10**7)))
cmax=30
cmin=30
small_start=$((-1*(10**5)))
large_start=$((-1*(10**8)))
final_time=$((2*(10**6))) 

for ((i=1; i<=8; i++)); do
   seed_array[i]=$((3696*i+2**i))
done

for ((i=2; i<=2; i++)); do
  for ((c=8; c<=8; c++)); do
  	./run.bash ${seed_array[$i]} \
  	"$run_name"_cycles"$c"_seed"${seed_array[$i]}"_AbsSin_"$Cmax"_"$Cmin"\
  	$vbar $dbar $c "FALSE" $Cmax $Cmin $small_time $large_time $final_time
  	
  #	./run.bash ${seed_array[$i]} \
  #	"$run_name"_cycles"$c"_seed"${seed_array[$i]}"_FullSin_"$Cmax"_"$Cmin"\
  #	$vbar $dbar $c "TRUE" $Cmax $Cmin $small_time $large_time $final_time  	
	done
done


