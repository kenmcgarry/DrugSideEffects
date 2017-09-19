# reviewers_JBI_validate_score.R
# validate the model created by using randomization and 10-fold cross-validation.
# score the models proposed camdidate drugs with respect to conventional drug parameters.

library(igraph)
library(dplyr)
library(Matrix)
library(tidyverse)
library(stringr)
library(ggplot2)
library(xtable)
library(ReactomePA)
library(clusterProfiler)
library(GOSemSim)
library(DOSE)


setwd("C:/R-files/sider")    # point to where my code lives
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData")
load("8thJune2017.RData")
source("reviewers_JBI_functions.R")  # 

# Jaccard similarity coefficient is used to generate an association score to calculate the
# overlap between two sets of features i.e. :
#   1. off-targets & pathway      2. drug chemical similarity & side-effects
#   3. on-targets & pathway       4. disease ontology pairs 
#   5. side-effects & drug ATC
#
# The Jaccard association index then draws together the combined individual scores for each pair and
# providing a value between 0 and 1 which is used to rank each candidate drug.
  
# Calculate on-target promiscuity score for all drugs
on_count <- length(unique(candidate_targets[,2])) # how many ontarget proteins in total?
on_score <- vector()
for (i in 1:nrow(candidate_list)){
  on_temp <- filter(candidate_targets, name == candidate_list[i,2])
  on_score[i] <- length(unique(on_temp$SYMBOL))
# how many ontargets attached to each drug?
}

on_score <- on_score/on_count   # ratio of on-targets to all on-targets for each candidate drug
rm(on_temp,on_count) # tidy up mess

# Calculate side-effect similarity score for all drugs
se_score <- mycandidates[,3]/100  # divide by 100 to remove percentage and get a value between 0 and 1


# Calculate score for number of diseases/symptoms each candidate drug treats
# inds_score
inds_score <- vector()
for (i in 1:length(inds)){
  inds_temp <- unlist(inds[i])
  inds_score[i] <- length(unique(inds_temp))
}
inds_count <- length(unique(unlist(inds)))
inds_score <- inds_score/inds_count


collate_scores <- cbind(inds_score,se_score,on_score)
collate_scores <- jaccard(collate_scores)  # ensure Matrix library is loaded.
collate_scores <- diag(collate_scores)
candidate_list_scores <- cbind(candidate_list,collate_scores)

head(candidate_list_scores)
candidate_list_scores <- arrange(candidate_list_scores,desc(collate_scores))
head(candidate_list_scores)

print.xtable(xtable(candidate_list_scores))


#==================================================================
# https://stats.stackexchange.com/questions/104040/resampling-simulation-methods-monte-carlo-bootstrapping-jackknifing-cross
# https://github.com/jrnold/resamplr
# https://drsimonj.svbtle.com/easy-leave-one-out-cross-validation-with-pipelearner

Yvar <- c(8,9,10,13,12, 14,18,12,8,9,   1,3,2,3,4)
Xvar <- c(rep("A", 5),  rep("B", 5),    rep("C", 5))
mydf <- data.frame (Yvar, Xvar)

boot.samples <- list()
for(i in 1:10) {
  t.xvar <- Xvar[ sample(length(Xvar), length(Xvar), replace=TRUE) ]
  t.yvar <- Yvar[ sample(length(Yvar), length(Yvar), replace=TRUE) ]
  b.df <- data.frame (t.xvar, t.yvar) 
  boot.samples[[i]] <- b.df 
}
str(boot.samples)
boot.samples[1]

permt.samples <- list()
for(i in 1:10) {
  t.xvar <- Xvar[ sample(length(Xvar), length(Xvar), replace=FALSE) ]
  t.yvar <- Yvar[ sample(length(Yvar), length(Yvar), replace=FALSE) ]
  b.df <- data.frame (t.xvar, t.yvar) 
  permt.samples[[i]] <- b.df 
}
str(permt.samples)
permt.samples[1]

# ============================

# https://www.r-bloggers.com/matrix-factorization/
require(recommenderlab)
require(sanealytics)
data(MovieLense) 

unroll_Vecs <- function (params, Y, R, num_users, num_movies, num_features) {
  # Unrolls vector into X and Theta
  # Also calculates difference between preduction and actual 
  
  endIdx <- num_movies * num_features
  
  X     <- matrix(params[1:endIdx], nrow = num_movies, ncol = num_features)
  Theta <- matrix(params[(endIdx + 1): (endIdx + (num_users * num_features))], 
                  nrow = num_users, ncol = num_features)
  
  Y_dash     <-   (((X %*% t(Theta)) - Y) * R) # Prediction error
  
  return(list(X = X, Theta = Theta, Y_dash = Y_dash))
}

J_cost <-  function(params, Y, R, num_users, num_movies, num_features, lambda, alpha) {
  # Calculates the cost
  
  unrolled <- unroll_Vecs(params, Y, R, num_users, num_movies, num_features)
  X <- unrolled$X
  Theta <- unrolled$Theta
  Y_dash <- unrolled$Y_dash
  
  J <-  .5 * sum(   Y_dash ^2)  + lambda/2 * sum(Theta^2) + lambda/2 * sum(X^2)
  
  return (J)
}

grr <- function(params, Y, R, num_users, num_movies, num_features, lambda, alpha) {
  # Calculates the gradient step
  # Here lambda is the regularization parameter
  # Alpha is the step size
  
  unrolled <- unroll_Vecs(params, Y, R, num_users, num_movies, num_features)
  X <- unrolled$X
  Theta <- unrolled$Theta
  Y_dash <- unrolled$Y_dash
  
  X_grad     <- ((   Y_dash  %*% Theta) + lambda * X     )
  Theta_grad <- (( t(Y_dash) %*% X)     + lambda * Theta )
  
  grad = c(X_grad, Theta_grad)
  return(grad)
}

# Now that everything is set up, call optim
print(
  res <- optim(par = c(runif(num_users * num_features), runif(num_movies * num_features)), # Random starting parameters
               fn = J_cost, gr = grr, 
               Y=Y, R=R, 
               num_users=num_users, num_movies=num_movies,num_features=num_features, 
               lambda=lambda, alpha = alpha, 
               method = "L-BFGS-B", control=list(maxit=maxit, trace=1))
)


scheme <- evaluationScheme(MovieLense, method = "split", train = .9,k = 1, given = 10, goodRating = 4) 

recommenderRegistry$set_entry(method="RSVD", dataType = "realRatingMatrix", fun=REAL_RSVD,
  description="Recommender based on Low Rank Matrix Factorization (real data).")













