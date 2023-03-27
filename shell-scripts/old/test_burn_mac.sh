# !/bin/zsh
# model run 

conda activate stan


basedir="./full-model/fire-sims/burns/${burn_mod}/"
cd ${basedir}
sttime=$(date +"%d%b%Y_%H%M")
model="stan/${burn_mod}_${params}"
outfile="csv-fits/${burn_mod}_${suffix}_${params}_${sttime}.csv"

# compile model and link c++ 
cmdstan_model ${model}
# run model with 3 chains
./${model} sample num_samples=1500 num_warmup=1500 num_chains=3 \
                  adapt init_buffer=300 term_buffer=200 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outfile} refresh=50 \
                  num_threads=3




method = sample (Default)
  sample
    num_samples = 1000 (Default)
    num_warmup = 1000 (Default)
    save_warmup = 0 (Default)
    thin = 1 (Default)
    adapt
      engaged = 1 (Default)
      gamma = 0.050000000000000003 (Default)
      delta = 0.80000000000000004 (Default)
      kappa = 0.75 (Default)
      t0 = 10 (Default)
      init_buffer = 75 (Default)
      term_buffer = 50 (Default)
      window = 25 (Default)
    algorithm = hmc (Default)
      hmc
        engine = nuts (Default)
          nuts
            max_depth = 10 (Default)
        metric = diag_e (Default)
        metric_file =  (Default)
        stepsize = 1 (Default)
        stepsize_jitter = 0 (Default)
    num_chains = 1 (Default)
id = 1 (Default)
data
  file =  (Default)
init = 2 (Default)
random
  seed = 325891136 (Default)
output
  file = output.csv (Default)
  diagnostic_file =  (Default)
  refresh = 100 (Default)
  sig_figs = -1 (Default)
  profile_file = profile.csv (Default)
num_threads = 1 (Default)