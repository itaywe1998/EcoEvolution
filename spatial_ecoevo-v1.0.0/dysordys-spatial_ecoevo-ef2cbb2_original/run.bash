#!/bin/bash 
R CMD BATCH ecoevo_main.R
grep "Error" ecoevo_main.Rout > err.log
