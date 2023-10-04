Three sub-directories :

1.  ecoevo : Stable version of modified code to run simple time gaps (\~10k years) with prospering evolution , short processes which the species survive
2.  modif : advanced modified versions

-   backup_beforeSpeciation_WithOrigins : data structure was shifted in order to include origins per cell in tabular data. rhs_eval and main source files only
-   backup_SpeciationWorks : Speciation application fully implied as an independent identity key for each species in each cell. Bound from above by max_species numbeer, see main notes.rhs_eval and main source files only Rest of directory contains essential support files like other major directories, currently (25.09.23) same as backup_SpeciationWorks

3.  original : Contains two main files.

-   ecoevo_main_original - The clone of Nornberg's code - runs short times for multiple species, for \~3 minutes, in contrast to most of other, shorter modified simulations

-   ecoevo_main - Current (25.09.23) working main, modified in attemp to find suitable separation between extinction on short adaptation time, and survival on long adaptation time.\

4.  (temporary) TempratureDebug :special directory with much resemblence to "original"'s file structure, created in order to detect dependencies in average temperature.

TODO : Decide if each important version should reside on its own branch on git, renaming directories to shorter and more clear names, and make order withing inner versions. Current sub-directories keep the separation over supporting files (all except rhs_eval and main, for instance plotting_functions stays the same in each major sub-directory but can differ between them.
