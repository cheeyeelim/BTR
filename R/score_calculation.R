#' @title Inner function of calculating Boolean model score
#' 
#' @description
#' This function calculates a final model score from pairwise and penalty scores.
#' 
#' @param x matrix vector. Pairwise scores computed by dist_measure().
#' @param bmodel S4 BoolModel object. Model to be evaluated.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' @param detail logical. Whether to give more details in score calculation. Default to FALSE.
m_score = function(x, bmodel, max_varperrule, detail=F)
{ 
  #(1) Calculate score for distance between states
  y = mean(apply(x, 2, min)) #best
  #y = mean(apply(x, 1, min))

  #(2) Calculate penalty term
  #(A) To penalise having too low or too high number of model states, when compared to the number of data states.
  #Ideally, the number of model states >= the number of data states.
  #abs(number of model states - number of data states)
  za = abs(ncol(x) - nrow(x)) / (nrow(x) * length(bmodel@target)) #best
  #za = abs(ncol(x) - nrow(x)) / (nrow(x))

  #(B) To penalise having too many variables in the rules.
  var_len = list() #combine act and inh rules for each variable.
  for(i in 1:length(bmodel@target))
  {
    var_len = c(var_len, list(c(bmodel@rule_act[[i]], bmodel@rule_inh[[i]])))
  }
  
  #calculate number and fraction of each variable.
  var_len = sapply(var_len, function(x) x[!grepl('0', x)]) #to remove '0'
  num_var = sapply(var_len, function(x) length(strsplit(paste(x, collapse=''), 'v')[[1]])-1 ) #to count the number of v[0-9]s terms. e.g. v1s&v2s will be counted as 2 terms.
  num_var[num_var<0] = 0 #if the rule for the variable is completely empty, has to set the negative values to 0.
  
  frac_var = (num_var-max_varperrule)
  zb_ind = num_var > max_varperrule
  zb = sum(frac_var[zb_ind])
    
  #(3) Calculate final score.
  #Specify the constants for each penalty term.
  a = 1
  b = 1

  f = y + a*za + b*zb #f ranges from 0 to infinity.
  
  if(detail)
  {
    output = c(f, y, za, zb)
    names(output) = c('f', 'y', 'za', 'zb')
  } else
  {
    output = c(f)
    names(output) = c('f')
  }
  return(output)
}

#' @title Calculating Boolean model score wrt to a dataset
#' 
#' @description
#' This function calculates a score for a Boolean model wrt to a dataset.
#' 
#' @param bmodel S4 BoolModel object. Model to be evaluated.
#' @param istate data frame. Must have only 1 row, which represents 1 initial state.
#' @param fcdata matrix. Represents the expression data df.
#' @param overlap_gene character vector. Specify which genes are present in both model and data inputs.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' @param detail logical. Whether to give more details in score calculation. Default to FALSE.
#' 
#' @export
calc_mscore = function(bmodel, istate, fcdata, overlap_gene, max_varperrule, detail=F)
{ 
  #(1) Simulate each of these models.
  mdata = simulate_model(bmodel, istate)

  #(2) Perform gene filtering on model state space.
  fmdata = filter_dflist(mdata, overlap_gene)
  fmdata = fmdata+0 #convert logical matrix into numeric matrix.
  
  if(class(fcdata)!='matrix')
  {
    fcdata = as.matrix(fcdata)
  }
  
  #(3) Score each model state wrt to data state.
  score_matrix = man_dist(fcdata, fmdata) #The first must be the data state. This returns a matrix of row=data, col=model.
  
  #(4) Calculate score for distance between states
  y = mean(apply(score_matrix, 2, min)) #best
  
  #(5) Calculate penalty term
  #(A) To penalise having too low or too high number of model states, when compared to the number of data states.
  #Ideally, the number of model states >= the number of data states.
  #abs(number of model states - number of data states)
  za = abs(ncol(score_matrix) - nrow(score_matrix)) / (nrow(score_matrix) * length(bmodel@target)) #best
  
  #(B) To penalise having too many variables in the rules.
  var_len = list() #combine act and inh rules for each variable.
  for(i in 1:length(bmodel@target))
  {
    var_len = c(var_len, list(c(bmodel@rule_act[[i]], bmodel@rule_inh[[i]])))
  }
  
  #calculate number and fraction of each variable.
  var_len = sapply(var_len, function(x) x[!grepl('0', x)]) #to remove '0'
  num_var = sapply(var_len, function(x) length(strsplit(paste(x, collapse=''), 'v')[[1]])-1 ) #to count the number of v[0-9]s terms. e.g. v1s&v2s will be counted as 2 terms.
  num_var[num_var<0] = 0 #if the rule for the variable is completely empty, has to set the negative values to 0.
  
  frac_var = (num_var-max_varperrule)
  zb_ind = num_var > max_varperrule
  zb = sum(frac_var[zb_ind])
  
  #(6 Calculate final score.
  #Specify the constants for each penalty term.
  a = 1
  b = 1
  
  f = y + a*za + b*zb #f ranges from 0 to infinity.
  
  if(detail)
  {
    output = c(f, y, za, zb)
    names(output) = c('f', 'y', 'za', 'zb')
  } else
  {
    output = c(f)
    names(output) = c('f')
  }
  
  return(output)
}

#' @title Calculates pairwise Manhattan distances between two matrices
#' 
#' @description
#' This function calculates pairwise Manhattan distances between two matrices.
#' 
#' @param x matrix
#' @param y matrix
man_dist = function(x, y)
{
  z = matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (k in 1:nrow(y)) {
    z[, k] <- colSums(abs(t(x) - y[k, ]))
  }
  return(z)
}

#' @title Calculate true positive, true negative, false positive and false negative
#' 
#' @description
#' This function calculates the true positive, true negative, false positive and false negative values from the adjacency matrices.
#' 
#' @param inf_mat matrix. It should be adjacency matrix of inferred network.
#' @param true_mat matrix. It should be adjacency matrix of true network.
#' 
#' @export
validate_adjmat = function(inf_mat, true_mat)
{
  output = rcpp_validate(inf_mat, true_mat)
  names(output) = c('tp', 'tn', 'fp', 'fn')
  
  return(output)
}

#' @title Calculate precision, recall, f-score, accuracy and specificity 
#' 
#' @description
#' This function calculates the precision, recall, f-score, accuracy and specificity from the output of validate_adjmat().
#' 
#' @param x integer vector. Vector output by validate_adjmat().
#' 
#' @export
calc_roc = function(x)
{
  #(1) Calculate the precision/positive predictive value.
  #probability that the disease is present when the test is positive
  #range from 0 to 1.
  p = unname(x['tp'] / (x['tp'] + x['fp']))
  
  #(2) Calculate the negative predictive value.
  #probability that the disease is not present when the test is negative
  #range from 0 to 1.
  np = unname(x['tn'] / (x['tn'] + x['fn']))
  
  #(3) Calculate the recall/sensitivity/true positive rate.
  #probability that a test result will be positive when the disease is present.
  #range from 0 to 1.
  r = unname(x['tp'] / (x['tp'] + x['fn']))
  
  #(4) Calculate the accuracy.
  a = unname((x['tp'] + x['tn']) / (x['tp'] + x['tn'] + x['fp'] + x['fn']))
  
  #(5) Calculate the specificity/true negative rate.
  #probability that a test result will be negative when the disease is not present.
  #range from 0 to 1.
  s = unname(x['tn'] / (x['tn'] + x['fp']))
  
  #(6) Calculate the positive likelihood ratio.
  #ratio between the probability of a positive test result given the presence of the disease and the probability of a positive test result given the absence of the disease.
  plr = r / (1 - s)
  
  #(7) Calculate the negative likelihood ratio.
  #ratio between the probability of a negative test result given the presence of the disease and the probability of a negative test result given the absence of the disease
  nlr = (1 - r) / s
  
  #(8) Calculate the f-score.
  #harmonic average of precision and recall.
  f = 2*p*r / (r+p)
  
  output = c(p, np, r, a, s, plr, nlr, f)
  names(output) = c('p', 'np', 'r', 'a', 's', 'plr', 'nlr', 'f')
  
  return(output)
}
