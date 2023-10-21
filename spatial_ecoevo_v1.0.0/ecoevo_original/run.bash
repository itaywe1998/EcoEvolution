#!/bin/bash

seed=$1
run_name=$2
y=$3
x=$4
cycles=$5
updown=$6
Cmax=$7 
Cmin=$8


model="Tdep"
small="TRUE"
dual="TRUE"
periodic="TRUE"
do_run=true
if $do_run
then
echo "Short Adaptation Time"
Rscript ecoevo_main.R $model $small $seed $run_name $y $x $cycles $updown $Cmax $Cmin
small="FALSE"
echo "Long Adaptation Time"
Rscript ecoevo_main.R $model $small $seed $run_name $y $x $cycles $updown $Cmax $Cmin
fi

echo "Demo"
Rscript demo.R $dual $run_name $y $x $cycles $periodic
