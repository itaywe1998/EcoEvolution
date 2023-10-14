#!/bin/bash

seed=$1
run_name=$2
y=$3
x=$4

model="Tdep"
small="TRUE"
dual="TRUE"
do_run=true
if $do_run
then
Rscript ecoevo_main.R $model $small $seed $run_name $y $x
small="FALSE"
Rscript ecoevo_main.R $model $small $seed $run_name $y $x
fi

Rscript demo.R $dual $run_name $y $x
