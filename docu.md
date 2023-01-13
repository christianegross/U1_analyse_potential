# Documentation
## Operation
### Determination of potential, renormalized anisotropy and $r_0$

The observables have to be generated with the [su2-package](https://github.com/urbach/su2), the analysis is based on the [hadron-package](https://github.com/HISKP-LQCD/hadron).
The program has to be run in the folder where the results for the observables are.
Run the analysis-files with the option --help to get an overview of all possible parameters that can be changed.
The scripts will produce plots, results in a csv-table and an RData-object that saves more data that is needed for further analysis.

### Determination of continuum limit

This is done with the predictbeta scripts, the results will be csv-files and some plots, as well as RDS data for the fitresults.

### TODO
- [x] Document remaining code
- [ ] make function to read nsave in
- [ ] make function to read in bootsamples
- [ ] clean up plot functions
- [ ] command line options for predictbeta, predictbetanaive