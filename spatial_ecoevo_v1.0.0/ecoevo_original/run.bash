#!/bin/bash

seed=$1
run_name=$2
vbar=$3
dbar=$4
cycles=$5
updown=$6
Cmax=$7 
Cmin=$8
small_time=$9
large_time=$10
final_time=$11


model="Tdep"
small="TRUE"
dual="TRUE"
periodic="TRUE"
do_run=true
if $do_run
then
echo "Short Adaptation Time"
Rscript ecoevo_main.R $model $small $seed $run_name $vbar $dbar $cycles 
                      $updown $Cmax $Cmin $small_time $final_time
small="FALSE"
echo "Long Adaptation Time"
Rscript ecoevo_main.R $model $small $seed $run_name $vbar $dbar $cycles 
                      $updown $Cmax $Cmin $large_time $final_time
fi

echo "Demo"
Rscript demo.R $dual $run_name $vbar $dbar $cycles $periodic 
               $small_time $large_time $final_time
