#!/bin/bash
IN=$1
arrIN=(${IN//#/ })

seed=${arrIN[0]}
run_name=${arrIN[1]}
vbar=${arrIN[2]}
dbar=${arrIN[3]}
cycles=${arrIN[4]}
updown=${arrIN[5]}
Cmax=${arrIN[6]} 
Cmin=${arrIN[7]}
small_time=${arrIN[8]}
large_time=${arrIN[9]}
final_time=${arrIN[10]}


model="Tdep"
small="TRUE"
dual="TRUE"
do_run=true

if $do_run
then
  echo "Short Adaptation Time"
  input="${model}#${small}#${seed}#${run_name}#${vbar}#${dbar}#${cycles}#${updown}#${Cmax}#${Cmin}#${small_time}#${final_time}"
  Rscript ecoevo_main.R $input
  small="FALSE"
  echo "Long Adaptation Time"
  input="${model}#${small}#${seed}#${run_name}#${vbar}#${dbar}#${cycles}#${updown}#${Cmax}#${Cmin}#${large_time}#${final_time}"
  Rscript ecoevo_main.R $input
fi

echo "Demo"
Rscript demo.R $dual $run_name $vbar $dbar $cycles
