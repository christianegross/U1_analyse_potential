# Workflow

This document lists the steps taken to generate input files and commands for taking a continuum limit in temporal direction as described in arXiv:2503.11480.

## Install libraries

Install the su2 library according to the installation instructions in the folder how_to_install.
We provide instructions for debian based machines and the typical modular structure of HPC-clusters.

The library provides functionalities for different gauge groups, here we will only consider the U(1) group.

## Generating configurations

Decide what configurations you want to generate and if you want to generate them with an mcmc or a heatbath algorithm.
Put the informations for the configs into files like paramsinput_contlimitL16.csv or paramsinput_contlimitL3.csv for the mcmc simulations and in files like paramsinput_heatbath.csv for the heatbath algorithm.

Run the scripts print_input_heatbath.sh or print_input_mcmc.sh top generate input files, make sure to adapt the variables for file storage to your system. 
Adapt confdirfile, resdirfile, inputfile and inputfileqbig.

Give the following informations:
- beta: coupling constant of the config
- Ns: spatial lattice sites
- Nt: temporal lattice sites
- xi_in: input anisotropy
- nape: number of APE smears, recommended 0
- alpha: alpha-parameter for APE smearing
- cores: number of cores to simulate the config
- betaone: put this only if you want to keep several simulations separate
- fraction: maximum extent of the Wilson loops that are measured, recommended 0.5
- skip: number of saved configurations to skip for thermalization
- meas: number of configurations that are simulated
- nsave: every nsaveth configurations is saved. Wilson loops are only measured on every nsaveth configueration
- offset: set this to 1 if you measure the wilson loops also at W(x=0, y=t, t=0)
- every: only use every "everyth" config for the analysis
- startheat: 0 for cold start, 1 for hot start
- n_overrelax: number of overrelaxation steps to take after every sweeps (only heatbath algorithm)
- n_heatbath: number of heatbath steps taken after every sweep (only heatbath algorithm)

## Analysing configurations

Run analysissubtracted for the normal potential and analysisrotated for the sideways potential with the appropriate parameters for every configuration.
When generating the inputfiles, the correct parameters are output automatically into the files commandsRmcmc.txt and commandsRheatbath.txt.

To draw the bootstrapsamples and compute the effective masses, first run the scripts with the options --analyse --drawbootstrap. 
Rerun with only --analyse to get more condensed plots of the efefctive mass fits.

Run with the option --dofit to determine the anisotropy and $r_0$ for the different analysis chains.

The results are written in the folder where the scripts are run, with a summary written in a resultsummary.csv file

## Determining continuum limit for the single analysis chains

If necessary, copy the different resultsummaryfiles together (one for normal, one for sideways), and copy or link the result files to the folder where you want to run the further analysis.

Run predictbeta.R, for example with the commands provided in contlimitcommands.sh, to produce an overview and plots to determine the matching beta.

Select the matching beta for each analysis chain and put the results in a table like inputpredictcontlim.csv.

Determine the continuum limit in the large colume by calling chosepredict.R, for example with the commands provided in contlimitcommands.sh.

Generate configurations in the small volume as in the preceding paragraph.
Analyse the data with L3singleplaquette.R. 
The commands are written to the commandsR files by the code that produces the input files.

Determine the continuum limit at small volume by calling L3contlimit.R, for example with the commands provided in contlimitcommands.sh.


 ## Determining combined continuum limit and matching level

 Run average_contlimits.r with the correct options, for example with the commands provided in contlimitcommands.sh.


 Run makescatterplotsfinalresult.R. 
 Adapt the locations of the files used in the script to your machine.

 Run collectfinalresult.Rmd.
 Adapt the locations of the files used in the script to your machine.
