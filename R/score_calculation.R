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
#' @param model_encoding numeric. Only required if the bmodel is encoded.
#' @param detail logical. Whether to give more details in score calculation. Default to FALSE.
#' 
#' @export
calc_mscore = function(bmodel, istate, fcdata, overlap_gene, max_varperrule, model_encoding, detail=F)
{ 
  if(class(bmodel)!='BoolModel')
  {
    bmodel = decompress_bmodel(bmodel, model_encoding)
  }
  
  #(1) Simulate each of these models.
  fmdata = simulate_model(bmodel, istate)
 
  #(2) Perform gene filtering on model state space.
  fmdata = filter_dflist(fmdata, overlap_gene)
  fmdata = fmdata+0 #convert logical matrix into numeric matrix.
  
  if(class(fcdata)!='matrix')
  {
    fcdata = as.matrix(fcdata)
  }

  #(3) Score each model state wrt to data state.
  score_matrix = man_dist(fcdata, fmdata) #The first must be the data state. This returns a matrix of row=data, col=model.
  
  #(4) Get the orderings of data states wrt to model states.
  rank_matrix = apply(score_matrix, 2, rank)
  rownames(rank_matrix) = seq(1,nrow(rank_matrix))
  
  rank_matrix = apply(rank_matrix, 2, function(x) as.numeric(names(sort(x))))
  
  ans_row = rank_matrix[1,]
  for(i in 1:nrow(rank_matrix))
  {
    if(!anyDuplicated(ans_row)) #stop if there is no more duplicates.
    {
      break
    }
    
    if(i != nrow(rank_matrix))
    {
      ans_row[duplicated(ans_row)] = rank_matrix[i+1,][duplicated(ans_row)]
    }
  }
  ans_row[duplicated(ans_row)] = rank_matrix[1,][duplicated(ans_row)] #if any duplicate remains at the end of loop, set them back to the best value.
  
  #(5) Calculate score for distance between states
  y = mean(sapply(1:length(ans_row), function(x) score_matrix[ans_row[x],x])) / ncol(fmdata)

  #(6) Calculate penalty term
  #(A) To penalise the states. The proportions of 0s and 1s are good predictors.
  ybin = infotheo::discretize(fmdata, nbins=2) #force mininum bin to 2.
  ybin = unique(ybin)
  
  tmp_y = rle(sort(unlist(ybin, use.names = F)))$lengths
  names(tmp_y) = rle(sort(unlist(ybin, use.names = F)))$values
  
  tmp_x = c(0.5, 0.5)
  tmp_y = tmp_y / sum(tmp_y)
  
  za = exp(-entropy::chi2.empirical(tmp_x, tmp_y))
  
  #(B) To penalise having too many variables in the rules.
  var_len = list() #combine act and inh rules for each variable.
  for(i in 1:length(bmodel@target))
  {
    var_len = c(var_len, list(c(bmodel@rule_act[[i]], bmodel@rule_inh[[i]])))
  }
  
  #calculate number and fraction of each variable.
  var_len = sapply(var_len, function(x) x[!grepl('^0$', x)]) #to remove '0'
  num_var = sapply(var_len, function(x) length(strsplit(paste(x, collapse=''), 'v')[[1]])-1 ) #to count the number of v[0-9]s terms. e.g. v1s&v2s will be counted as 2 terms.
  num_var[num_var<0] = 0 #if the rule for the variable is completely empty, has to set the negative values to 0.
  
  frac_var = num_var-max_varperrule #only punished if above set threshold.
  frac_var[frac_var<0] = 0
  zb = mean(frac_var/max_varperrule) 
  
  #(7) Calculate final score.
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
    z[,k] = colSums(abs(t(x) - y[k,]))
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
