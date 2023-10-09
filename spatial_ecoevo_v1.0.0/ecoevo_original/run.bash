#!/bin/bash

model="Tdep"
small="TRUE"
seed=1234
run_name="DualDisp&Run"
y=200
dual="TRUE"
do_run=false

if $do_run
then
Rscript ecoevo_main.R $model $small $seed $run_name $y
small="FALSE"
Rscript ecoevo_main.R $model $small $seed $run_name $y
fi

Rscript demo.R $dual $run_name $y
