#!/bin/bash

seed_array=()
run_name="PeriodicCC_Factory"
vbar=$(bc -l <<<"3/(10^5)") #sort of scientific writing (not supported directly in bash)
dbar=$(bc -l <<<"1/(10^7)")
Cmax=25
Cmin=10
small_time=$((-1*(10**5)))
large_time=$((-1*(10**8)))
final_time=$((2*(10**6))) 

for ((i=1; i<=16; i++)); do
   seed_array[i]=$((3688*i+2**i))
done

for ((i=9; i<=16; i++)); do
  for ((c=4; c<=4; c++)); do
  	./run.bash "${seed_array[$i]}#"$run_name"_cycles"$c"_seed"${seed_array[$i]}"_AbsSin_"$Cmax"_"$Cmin"#$vbar#$dbar#$c#FALSE#$Cmax#$Cmin#$small_time#$large_time#$final_time"
	done
done


