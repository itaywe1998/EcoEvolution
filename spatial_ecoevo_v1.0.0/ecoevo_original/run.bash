#!/bin/bash
model="Tdep"
small="TRUE"
seed=1234
run_name="DualRunCheck"
Rscript ecoevo_main.R $model $small $seed $run_name
small="FALSE"
Rscript ecoevo_main.R $model $small $seed $run_name
