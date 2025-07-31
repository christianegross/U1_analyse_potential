# Analysis of the static potential

The scripts in this repository can be used to determine the static potential, the renormalized anisotropy $\xi_R$, the renormalized $\beta$ and the Hamiltonian limit of the plaquette

There is further information in the how_to_install folder and the example folder.

## Scripts

- myfunctions.R
- potential/ single ensemble:
    - analysisrotated.R
    - analysissubtracted.R
    - L3singleplaquette.R
- results continuum limits:
    - predictbeta.R
    - chosepredict.R (example inputpredictcontlim.csv)
    - L3contlimit.R (example inputpredictcontlim.csv)
- further analysis continuum limits:
    - average_contlimits.r
    - makescatterplotsfinalresults
    - collectfinalresult.Rmd
- helpers:
    - determine_correlation_contlim.R
    - matchwithellipse.R

There are various other scripts to deal with other analysis possibilities of the single ensembles, and other ways to ensure constant lattice spacing/take the continuum limit, and various scripts to format tables or do (intermediate) plots.

### Required packages
 - example installation scripts are in how_to_install: one script for debian based machines, tested on linux mint cinnamon 22.1, and one script for hpc clusters tested on juwels-booster at FZ Juelich and Marvin at Uni Bonn.
 - hadron and its dependencies -> statistical analysis, specifically needs commit dec9d57
 - optparse and its dependencies -> command line parameters

### Workflow

A more detailed workflow is in [the example folder](https://github.com/christianegross/U1_analyse_potential/tree/clean_for_publication/examples/README.md)

 - generate a collections of ensembles at the large volume with the [su2 code](https://github.com/HISKP-LQCD/su2). Select branch su3 or nested_sampling to enable heatbath simulations.
 - analyse them with analysisrotated/analysissubtracted
 - analyse those results with predictbeta
 - choose the matching beta or go back to the first step to generate more ensembles if needed
 - run the simulations at the matching beta at small volumes, analyse with L3singleplaquette
 - run chosepredict and after that L3contlimit
 - run the scripts for further analysis in the order in the list

