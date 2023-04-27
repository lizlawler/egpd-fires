data {
  int<lower=1> R; // # of regions
  int<lower=1> p; // # of parameters
  int<lower=1> T_all; // # of timepoints in full dataset
  int<lower=1> T_train;
  int<lower=1> T_hold;

  // covariate data
  array[R] matrix[T_all, p] X_full; // design matrix; 1-D array of size R with matrices T x p
  array[R] matrix[T_train, p] X_train; // design matrix; 1-D array of size R with matrices T x p
  
  // lower bound of burns
  real y_min;

  // training data
  int<lower=1> N_tb_obs;
  int<lower=1> N_tb_mis;
  int<lower=1> N_tb_all;
  array[N_tb_obs] real<lower=y_min> y_train_obs; // burn area for observed training timepoints
  array[N_tb_obs] int<lower=1> ii_tb_obs;
  array[N_tb_mis] int<lower=1, upper=N_tb_all> ii_tb_mis;
  array[N_tb_all] int<lower=1, upper=N_tb_all> ii_tb_all; // for broadcasting
  array[T_train] int<lower=1> idx_train_er;

  // holdout data
  int<lower=1> N_hold_obs;
  int<lower=1> N_hold_all; // includes 'missing' and observed
  array[N_hold_obs] int<lower=1> ii_hold_obs; // vector of indices for holdout data timepoints
  array[N_hold_all] int<lower=1> ii_hold_all; // vector of indices for broadcasting to entire holdout dataset
  array[N_hold_obs] real<lower=1> y_hold_obs; // burn area for observed holdout timepoints
  array[T_hold] int<lower=1> idx_hold_er;

  // neighbor information
  int<lower=0> n_edges;
  array[n_edges] int<lower=1, upper=R> node1; // node1[i] adjacent to node2[i]
  array[n_edges] int<lower=1, upper=R> node2; // and node1[i] < node2[i]

  // indicator matrices for ecoregions
  matrix[R, R] l3;
  matrix[R, R] l2;
  matrix[R, R] l1;

  // indicator matrices for AR(1) penalization on spline coefficients of betas
  matrix[p, p] equal;
  matrix[p, p] bp_lin;
  matrix[p, p] bp_square;
  matrix[p, p] bp_cube;
  matrix[p, p] bp_quart;
}