#!/bin/bash 
R CMD BATCH ecoevo_main.R
grep "Hello" ecoevo_main.Rout > Ancestry.log
grep "Error" ecoevo_main.Rout > err.log
