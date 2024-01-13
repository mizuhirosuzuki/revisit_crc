# Replication package for "Revisiting Correlated Random Coefficient Model in Technology Adoption"

This repo stores codes for reproducing the result in my paper "Revisiting Correlated Random Coefficient Model in Technology Adoption".

There are two R scripts in the repo (both in the `R` folder):
- `CRC_no_covariate.R`: a script to reproduce the results in the case where no covariate is included in the regressions;
- `CRC.R`: a script to reproduce the results in the case where a covariate is included in the regressions.

(I don't have an `renv` setup in this repo... sorry!)

To reproduce the results, you can use the `Makefile` of this repo:
- `make simulation`: run the two R scripts to reproduce the results and save them to `Figures` folder;
- `all`: run the two R scripts for the results and compile the `.tex` file to get a pdf of the paper.

Please feel free to contact me if you find there is any issue.
You can leave a message on the issue page or send an email to me (mizuhiro.suzuki@gmail.com).
