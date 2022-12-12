# Analysis of the static potential

The scripts in this repository can be used to determine the static potential, the renormalized anisotropy $\xi_R$, the renormalized $\beta$ and the Hamiltonian limit of the plaquette

## Scripts

- myfunctions.R
- potential:
    - analysissmall.R
    - analysisrotated.R
    - analysissubtracted.R
- results:
    - predictbeta.R
    - predictbetanaive.R
- helpers:
    - makeinputtable.R
    - converttablelatex.R

## Operation
### Determination of potential, renormalized anisotropy and $r_0$

The observables have to be generated with the [su2-package](https://github.com/urbach/su2), the analysis is based on the [hadron-package](https://github.com/HISKP-LQCD/hadron).
The program has to be run in the folder where the results for the observables are.
Run the analysis-files with the option --help to get an overview of all possible parameters that can be changed.
The scripts will produce plots, results in a csv-table and an RData-object that saves more data that is needed for further analysis.

### Determination of continuum limit

This is done with the predictbeta scripts, the results will be a csv-file and some plots