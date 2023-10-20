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
do_run=true
if $do_run
then
Rscript ecoevo_main.R $model $small $seed $run_name $y $x $cycles $updown $Cmax $Cmin
small="FALSE"
Rscript ecoevo_main.R $model $small $seed $run_name $y $x $cycles $updown $Cmax $Cmin
fi

Rscript demo.R $dual $run_name $y $x
