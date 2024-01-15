data {
  int<lower=1> R; // # of regions
  int<lower=1> p; // # of parameters
  int<lower=1> T_all; // # of timepoints in full dataset
  int<lower=1> T_train;
  int<lower=1> T_hold;
  
  // covariate data
  array[R] matrix[T_all, p] X_full; // design matrix; 1-D array of size r with matrices t x p
  array[R] matrix[T_train, p] X_train;
  
  // area offset
  array[R] real area_offset; // known offset vector of areas
  
  // training dataset
  array[T_train] int<lower=1> idx_train_er; // vector of indices for training data timepoints
  array[T_train, R] int<lower=0> y_train_count; // response data
  
  // holdout dataset
  array[T_hold] int<lower=1> idx_hold_er; // vector of indices for holdout data timepoints
  array[T_hold, R] int<lower=0> y_hold_count; // response data
  
  // neighbor information
  int<lower=0> n_edges;
  array[n_edges] int<lower=1, upper=R> node1; // node1[i] adjacent to node2[i]
  array[n_edges] int<lower=1, upper=R> node2; // and node1[i] < node2[i]
  
  // indicator matrices for ecoregions
  matrix[R, R] l3;
  matrix[R, R] l2;
  matrix[R, R] l1;
  
  // indicator matrices for AR(1) process on betas
  matrix[p, p] equal;
  matrix[p, p] bp_lin;
  matrix[p, p] bp_square;
  matrix[p, p] bp_cube;
  matrix[p, p] bp_quart;
}