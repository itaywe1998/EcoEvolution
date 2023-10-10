#!/bin/bash

seed=$1
run_name=$2
y=$3

model="Tdep"
small="TRUE"
dual="TRUE"
do_run=true
if $do_run
then
Rscript ecoevo_main.R $model $small $seed $run_name $y
small="FALSE"
Rscript ecoevo_main.R $model $small $seed $run_name $y
fi

Rscript demo.R $dual $run_name $y
