# EGPD - Wildfires in the contiguous United States

Code and data associated with Lawler and Shaby's "Anthropogenic and meteorological effects on the counts and sizes of moderate and extreme wildfires" (in review).

## Table of contents

-   [Data processing](#data-processing)
-   [Modeling wildfire counts and sizes](#modeling-wildfire-counts-and-sizes)
-   [Results and figures](#results)

## Data processing {#data-processing}

-   The R, python, and shell scripts found in [data/](/data/) are numbered and should be executed in chronological order. Any two scripts with the same number indicate that they do not rely on each other and can therefore be run in parallel.
-   These scripts begin with downloading the raw data and proceed through the processing, cleaning, and aggregation of the various datasets for use in the models.
-   Should you wish to skip the downloading, processing, and cleaning steps, you can use the included files in [`data/processed/`](data/processed/) to execute [`data/04_make_stan_data_lists.R`](data/04_make_stan_data_lists.R), which creates the .json and .RDS files found in [`data/stan_lists/`](data/stan_lists/). Alternatively, you can start straight from the the lists of data created for the Stan models.

## Modeling wildfire counts and sizes {#modeling-wildfire-counts-and-sizes}

All models were run using `cmdstan`, which is executed from the command line. Code to fully execute these models on an HPC with slurm scheduling exists in the [`shell_scripts/slurm_alpine`](shell_scripts/slurm_alpine/) folder. Code to run some of these models on a remote Linux machine exists in the [`mp_server`](shell_scripts/mp_server) folder.

### Burned sizes submodel: phased approach

Stan code has been created for every permutation of the EGPD carrier families and lognormal distribution for the burned sizes submodels. However, we ran these models in a phased approach given the sheer number of combinations when incorporating all permutations of the parameters and datasets.

#### Phase one

-   Execute [`shell_scripts/slurm_alpine/run_g1_params_NUTS.sh`](shell_scripts/slurm_alpine/run_g1_params_NUTS.sh).

-   Once the models in phase one have finished running, proceed with 'Phase one' in the [model scoring](scores_traceplots/model_comparison.R) script.

#### Phase two

-   Execute [`shell_scripts/slurm_alpine/run_g1_datasets_NUTS.sh`](shell_scripts/slurm_alpine/run_g1_datasets_NUTS.sh).

-   Once the models have finished running, proceed with 'Phase two' in the [model scoring](scores_traceplots/model_comparison.R) script.

#### Phase three

-   Execute [`run_g2_NUTS.sh`](shell_scripts/slurm_alpine/run_g2_NUTS.sh), [`run_g3_NUTS.sh`](shell_scripts/slurm_alpine/run_g3_NUTS.sh), [`run_g4_NUTS.sh`](shell_scripts/slurm_alpine/run_g4_NUTS.sh), and [`run_lognorm_NUTS.sh`](shell_scripts/slurm_alpine/run_lognorm_NUTS.sh). NB: G3 and G4 take longer than seven days to complete, so will likely be unable to finish on an HPC with slurm scheduling. Please look to [`mp_server`](shell_scripts/mp_server) for guidance on running these locally.

-   Once the models have finished running, proceed with 'Phase three' in the [model scoring](scores_traceplots/model_comparison.R) script.

### Occurrences submodel: phased approach

-   Execute [`shell_scripts/slurm_alpine/run_counts_models_NUTS.sh`](shell_scripts/slurm_alpine/run_counts_models_NUTS.sh).

-   Once the six models have finished running, proceed with 'Phase one' in the [model scoring](scores_traceplots/model_comparison.R) script.

-   Execute [`shell_scripts/slurm_alpine/run_counts_datasets_NUTS.sh`](shell_scripts/slurm_alpine/run_counts_datasets_NUTS.sh), then proceed with 'Phase two' of scoring.

### Joint model

-   Execute [`shell_scripts/slurm_alpine/run_joint_NUTS.sh`](shell_scripts/slurm_alpine/run_joint_NUTS.sh) , which runs the sampler and then runs [`figures/extract_joint_mcmc_draws.R`](figures/extract_joint_mcmc_draws.R) .

## Results and figures

All R scripts needed to recreate the figures presented in the paper are included in the [`figures/`](figures/) folder. The MCMC draws and model scores must be extracted prior to running these scripts. The scripts can then be run in no particular order.
