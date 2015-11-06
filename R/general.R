#' @title Verbose cat
#' 
#' @description
#' This function is a simple wrapper for cat with a Boolean value for turning it on/off.
#' 
#' @param string character vector. String intended to be printed.
#' @param bool logical. Specify whether to print the string to stout or not.+
vcat = function(string, bool)
{
  if(bool)
  {
    cat(string)
  }
  
  invisible()
}

#' @title Extract Boolean terms
#' 
#' @description
#' This function extracts the terms within a Boolean rule, using OR as the separator. Bracketed variables are counted as one term.
#' 
#' @param brule character vector. It should be either an activating or inhibiting rule for a target gene.
extract_term = function(brule)
{
  if(brule!='0') #account for cases with empty rule.
  {
    #Strip outermost brackets.
    brule = gsub('^[(](.+)[)]$', '\\1', brule)
    term = strsplit(brule, '[|]')[[1]]
    term = gsub('^[(](.+)[)]$', '\\1', term)
  } else
  {
    term = '0'
  }

  return(term)
}

#' @title Check for matching terms
#' 
#' @description
#' This function checks if the first term is found within a vector of terms. Return a logical value.
#' It is smarter than simple string matching, e.g. match_term(v1s&v2s, v2s&v1s) == T, match_term(v1s, v1s&v2s) == T, match_term(v1s&v2s, v1s) == T.
#' 
#' @param t1 character vector of length 1 or more. It should be a vector of gene variable.
#' @param t2 character vector of length 1 or more. It should be a vector of gene variables.
#' @param mode character. Indicates the mode of action. Options: 'logic', 'unique'. Default to 'logic'.
match_term = function(t1, t2, mode='logic')
{
  if(t1[1] != '0')
  {
    t1 = unlist(strsplit(t1, '&')) #this will split up terms with &.
    
    #Check if t1 is present in t2.
    #Note: if t1='v1s', this will match 'v1s', 'v2s&v1s'.
    if(length(t2) > 1)
    {
      check_mat = sapply(t1, function(x) grepl(x, t2))
    } else
    {
      check_mat = t(matrix(sapply(t1, function(x) grepl(x, t2)))) #matrix() is necessary when t2 has length 1.
    }
    
    if(mode == 'logic')
    {
      check_ind = which(rowSums(check_mat) > 0)
      
      if(length(check_ind) == 0)
      {
        match = FALSE
      } else
      {
        match = TRUE
      }
      
      return(match)
    } else if(mode == 'unique')
    {
      check_ind = which(rowSums(check_mat) == 0)
      
      return(t2[check_ind])
    } else
    {
      stop('Error in specifying match_term() mode.')
    }
  } else if(t1[1] == '0')
  {
    if(mode == 'logic')
    {
      return(FALSE)
    } else if(mode == 'unique')
    {
      return(t2)
    } else
    {
      stop('Error in specifying match_term() mode.')
    }
  }
}

#' @title Check for matching states
#' 
#' @description
#' This function finds a match between two df of states. Returns a row index vector indicating for each row of mstate, what is the corresponding row in xstate. If a match cannot be found, a 0 will be return.
#' Only columns that are present in both df will be used in comparison. Note that the row index starts from 1 (as in R), not from 0 (as in cpp).
#' 
#' @param mstate data frame. It should be a state(row) x gene(column) df. colnames will be used in comparison.
#' @param xstate data frame. It should be a state(row) x gene(column) df. colnames will be used in comparison.
match_state = function(mstate, xstate)
{
  #Filtering the columns in mstate and xstate.
  same_col = intersect(colnames(mstate), colnames(xstate))
  fmstate = mstate[, same_col]
  fxstate = xstate[, same_col]
  
  ind = match_state_loop(as.matrix(fmstate), as.matrix(fxstate))
  
  return(ind) #zeroes are present due to mismatching in cpp code.
}

#' @title Pick a random minimum value
#' 
#' @description
#' This function locates the minimum value in a vector (similar to which.min), however it will randomly break ties when there are multiple minimum values.
#' 
#' @param x numeric vector.
#' @param favour_first logical. If this is TRUE, and the first value in the vector supplies is one of the min, the index for the first value will always be returned.
which.random.min = function(x, favour_first=F)
{
  id = which(x == min(x)) #get indices of all minimum values.
  
  # If the minimum value appear multiple times pick one index at random, otherwise just return its position in the vector
  # Exception occurs if the initial model (at position 1) also has the minimum value, then the initial model is favoured.
  if(length(id) > 1 & 1 %in% id & favour_first==T)
  {
    id = 1
  } else if(length(id) > 1) #for all other cases with multiple minimum values, pick a random one.
  {
    id = sample(id, 1)
  } else if(length(id) == 1) #if there is only single min value, return it.
  {
    id = id
  }
  
  return(id)
}

#' @title Filter columns of df in a list
#' 
#' @description
#' This function filters columns of multiple df in a list, when compared using a vector. Use through lapply().
#' 
#' @param x data frame. It should be a list element if used from lapply.
#' @param y character vector. It should contains gene names found in the colnames of the data frame.
#' @param uniq_bool logical. Whether to return unique rows only.
filter_dflist = function(x, y, uniq_bool=T)
{
  stopifnot(class(x)== 'data.frame' | class(x)== 'matrix')
  
  if(uniq_bool)
  {
    return(unique(x[,colnames(x) %in% y, drop=F]))
  } else
  {
    return(x[,colnames(x) %in% y, drop=F])
  }
}

#' @title Check if the Boolean model violates constraints.
#' 
#' @description
#' This function checks if the Boolean model violates contraints. Return logical value.
#' (1) Each gene rule should not have more terms than max_varperrule.
#' (2) The same term should not occur twice in the same rule.
#' 
#' @param bmodel S4 BoolModel object.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' 
#' @export
check_bmodel = function(bmodel, max_varperrule)
{
  check_bool = F
  
  act_rule = lapply(bmodel@rule_act, function(x) unlist(strsplit(x, '&')))
  inh_rule = lapply(bmodel@rule_inh, function(x) unlist(strsplit(x, '&')))
  
  #(1) Check 1 : The same term should not occur twice in the same rule.
  check_bool = check_bool | any(sapply(act_rule, function(x) any(duplicated(x))))
  check_bool = check_bool | any(sapply(inh_rule, function(x) any(duplicated(x))))
  
  #(2) Check 2 : Each gene rule should not have more terms than max_varperrule.
  act_rule = sapply(act_rule, function(x) unique(x))
  act_rule = sapply(act_rule, function(x) x[!(x %in% '0')]) #remove zeroes.
  
  inh_rule = sapply(inh_rule, function(x) unique(x))
  inh_rule = sapply(inh_rule, function(x) x[!(x %in% '0')]) #remove zeroes.
  
  check_bool = check_bool | any(sapply(act_rule, function(x) length(x)>max_varperrule))
  check_bool = check_bool | any(sapply(inh_rule, function(x) length(x)>max_varperrule))
  
  if(check_bool)
  {
    return(FALSE)
  } else
  {
    return(TRUE)
  }
}

#' @title Obtain parameters for bimodal distribution from real data
#' 
#' @description
#' This function obtains parameters for bimodal distribution. Returns 4 parameters: mu1, mu2, sig1, sig2.
#' 
#' @param x matrix. Input expression data. Col-genes, row-samples.
#' @param data_type character. Specify data types: qpcr, rnaseq.
#' 
#' @export
param_bimodal = function(x, data_type='qpcr')
{
  require(MASS)
  
  #(1) Initialise data.
  tmp = initialise_raw_data(x, data_type)
  x_con = tmp[[1]] #continuous
  x_bin = tmp[[2]] #discrete
  
  #rescale the data to remove zeroes.
  x_con = x_con + 0.0001
  
  all_mu1 = c()
  all_mu2 = c()
  all_sig1 = c()
  all_sig2 = c()
  for(i in 1:ncol(x_bin))
  {
    #(2) Extract the parameters for the two modal distributions. 
    #First modal - low expression, 0s. Second modal - high expression, 1s.
    x_lowmode = x_con[x_bin[,i]!=1, i]
    x_highmode = x_con[x_bin[,i]==1, i]

    #(3) Estimate parameters
    param1 = MASS::fitdistr(x_lowmode, 'lognormal')
    param2 = MASS::fitdistr(x_highmode, 'lognormal')
    
    #For checking.
    #hist(rlnorm(1000, param1$estimate[1], param1$estimate[2]))
    #hist(rlnorm(1000, param2$estimate[1], param2$estimate[2]))
    
    all_mu1 = c(all_mu1, param1$estimate[1])
    all_mu2 = c(all_mu2, param2$estimate[1])
    all_sig1 = c(all_sig1, param1$estimate[2])
    all_sig2 = c(all_sig2, param2$estimate[2])
  }
  
  return(list(all_mu1=all_mu1, all_mu2=all_mu2, all_sig1=all_sig1, all_sig2=all_sig2))
}


#' @title Generate random real numbers from binary values
#' 
#' @description
#' This function generates random real numbers from binary values, with supplied parameters. Returns a vector of real values.
#' 
#' @param x logical or 0/1 numeric matrix. Col-genes, row-samples.
#' @param param list of parameters given by param_bimodal().
#' 
#' @export
bin_to_real = function(x, param)
{ 
  require(MASS)
  
  #(1) Convert logical to numeric.
  x = x + 0
  
  #(2) Estimate the distribution for the parameters.
  mu1_dist = MASS::fitdistr(-param$all_mu1, 'lognormal')
  mu2_dist = MASS::fitdistr(-param$all_mu2, 'lognormal')
  sig1_dist = MASS::fitdistr(param$all_sig1, 'lognormal')
  sig2_dist = MASS::fitdistr(param$all_sig2, 'lognormal')
  
  y = matrix(NA, ncol=ncol(x), nrow=nrow(x))
  for(i in 1:ncol(x))
  {
    #(3) Generating random values from the distribution.
    mu1_est = -rlnorm(1, mu1_dist$estimate[1], mu1_dist$estimate[2])
    mu2_est = -rlnorm(1, mu2_dist$estimate[1], mu2_dist$estimate[2])
    sig1_est = rlnorm(1, sig1_dist$estimate[1], sig1_dist$estimate[2])
    sig2_est = rlnorm(1, sig2_dist$estimate[1], sig2_dist$estimate[2])
    
    #(4) For each gene, generate random expression values using the obtained random parameters.
    for(j in 1:nrow(x))
    {
      if(x[j,i]==0)
      {
        y[j,i] = rlnorm(1, mu1_est, sig1_est)
        
        if(y[j,i] > 1)
        {
          y[j,i] = 1
        } else if(y[j,i] < 0)
        {
          stop('Error in generating continuous values.')
        }
      } else
      {
        y[j,i] = rlnorm(1, mu2_est, sig2_est)
        
        if(y[j,i] > 1)
        {
          y[j,i] = 1
        } else if(y[j,i] < 0)
        {
          stop('Error in generating continuous values.')
        }
      }
    }
  }
  
  return(y)
}

#' @title Check for equivalent models
#' 
#' @description
#' This function checks if the two models have the same rules. Return a logical value. Only TRUE if each rule for each gene is the same.
#' 
#' @param bmodel1 S4 BoolModel object.
#' @param bmodel2 S4 BoolModel object.
#' @param inter_bool logical. Indicate whether to consider AND terms.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' 
#' @export
equi_model = function(bmodel1, bmodel2, inter_bool, max_varperrule)
{
  stopifnot(length(bmodel1@target)==length(bmodel2@target))
  stopifnot(get_encodings(bmodel1, inter_bool)==get_encodings(bmodel2, inter_bool))
  
  ind = get_encodings(bmodel1, inter_bool)
  
  dist = model_dist(bmodel1, bmodel2, inter_bool, max_varperrule)
  
  if(dist==0)
  {
    match = T
  } else
  {
    match = F
  }
  
  return(match)
}

#' @title Calculate distance between Boolean models
#' 
#' @description
#' This method takes in two models and calculate the distance between them. The value return indicate the number of steps between the two models.
#' 
#' @param x S4 BoolModel object. Test model.
#' @param y S4 BoolModel object. Reference model.
#' @param inter_bool logical. Indicate whether to consider AND terms.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' 
#' @export
model_dist = function(x, y, inter_bool, max_varperrule)
{
  set_diff = unlist(model_setdiff(x, y, inter_bool, max_varperrule))
  
  #Calculate total dist.
  t_dist = length(set_diff)
  
  #Return results.
  return(t_dist)
}

#' @title Find the set difference between two Boolean models
#' 
#' @description
#' This method takes in two models and find the set difference between them. Return a vector with the set difference.
#' 
#' @param x S4 BoolModel object. Test model.
#' @param y S4 BoolModel object. Reference model.
#' @param inter_bool logical. Indicate whether to consider AND terms.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' @param directed logical. If TRUE, return the difference in terms with respect to x.
#' 
#' @export
model_setdiff = function(x, y, inter_bool, max_varperrule, directed=F)
{
  stopifnot(length(x@target) == length(y@target))
  stopifnot(get_encodings(x, inter_bool)==get_encodings(y, inter_bool))
  
  ind = get_encodings(x, inter_bool)
  
  x1 = compress_bmodel(x, ind, max_varperrule)
  x2 = compress_bmodel(y, ind, max_varperrule)
  
  #Pick which model has more terms, and which has less.
  more_ind = which.max(c(length(x1), length(x2)))
  less_ind = ifelse(more_ind==1, 2, 1)
  
  more_x = list(x1, x2)[[more_ind]]
  less_x = list(x1, x2)[[less_ind]]
  
  #Remove empty terms.
  more_x = more_x[more_x%%1!=0]
  less_x = less_x[less_x%%1!=0]
  
  #Remove 0 and !0.
  #Find the encodings for 0 and !0.
  zero_enc = as.numeric(paste('0.', ind[which(names(ind)=='0')], sep=''))
  notzero_enc = as.numeric(paste('0.', ind[which(names(ind)=='!0')], sep=''))
  
  #Remove them.
  more_x = more_x[round(more_x-round(more_x), 4)!=zero_enc]
  more_x = more_x[round(more_x-round(more_x), 4)!=notzero_enc]
  
  less_x = less_x[round(less_x-round(less_x), 4)!=zero_enc]
  less_x = less_x[round(less_x-round(less_x), 4)!=notzero_enc]
  
  #Check for set difference.
  more_bool = sapply(more_x, function(x) !(x %in% less_x))
  less_bool = sapply(less_x, function(x) !(x %in% more_x))
  
  more_diff = more_x[more_bool]
  less_diff = less_x[less_bool]
  
  more_vec = c()
  if(length(more_diff)>0)
  {
    gene = round(more_diff)
    rule = gsub('^[0-9]+[.]', '', gsub(' ', '', format(more_diff, nsmall=4)))
    
    gene = x@target[gene]
    rule = unname(sapply(rule, function(x) names(which(ind==x))[1]))
    
    more_vec = rule
    names(more_vec) = gene
  }
  
  less_vec = c()
  if(length(less_diff)>0)
  {
    gene = round(less_diff)
    rule = gsub('^[0-9]+[.]', '', gsub(' ', '', format(less_diff, nsmall=4)))
    
    gene = x@target[gene]
    rule = unname(sapply(rule, function(x) names(which(ind==x))[1]))
    
    less_vec = rule
    names(less_vec) = gene
  }
  
  #Handle step difference that includes changing to and from AND terms.
  if(length(more_vec) != 0 & length(less_vec) != 0)
  {
    for(i in 1:length(more_vec))
    {
      for(j in 1:length(less_vec))
      {
        if(names(more_vec[i])==names(less_vec[j]))
        {
          if(grepl('!', more_vec[i]) & grepl('!', less_vec[j])) #for both inh terms. do not cross compare act and inh terms.
          { 
            #Split & terms if present.
            mvec = unlist(strsplit(more_vec[i], '&'))
            lvec = unlist(strsplit(less_vec[j], '&'))
            
            #Removes !.
            mvec = gsub('!', '', mvec)
            lvec = gsub('!', '', lvec)
            
            #Check for presence of &.
            mand_bool = ifelse(grepl('&', more_vec[i]), T, F)
            land_bool = ifelse(grepl('&', less_vec[j]), T, F)
            
            mbool = mvec %in% lvec
            lbool = lvec %in% mvec
            
            #assign NA if there is no term remaining.
            if(length(mvec[!mbool])!=0)
            {
              if(!is.na(mvec[!mbool][1]))
              {
                more_vec[i] = paste(mvec[!mbool], collapse='&')
              } 
            } else
            {
              more_vec[i] = NA
            }
            
            if(length(lvec[!lbool])!=0)
            {
              if(!is.na(lvec[!lbool][1]))
              {
                less_vec[j] = paste(lvec[!lbool], collapse='&')
              } 
            } else
            {
              less_vec[j] = NA
            }
            
            #if this involves AND term, adds & indicator to term name.
            if(mand_bool & !is.na(more_vec[i]) & !grepl('&', more_vec[i]))
            {
              more_vec[i] = paste('&', more_vec[i], sep='')
            }
            if(land_bool & !is.na(less_vec[j]) & !grepl('&', less_vec[j]))
            {
              less_vec[j] = paste('&', less_vec[j], sep='')
            }
            
            #Add the ! back to the terms.
            if(!is.na(more_vec[i]))
            {
              if(grepl('^&', more_vec[i]))
              {
                more_vec[i] = gsub('(^[&])(.+)', '\\1!\\2', more_vec[i])
              } else
              {
                more_vec[i] = paste('!', more_vec[i], sep='')
              }
            }
            if(!is.na(less_vec[j]))
            {
              if(grepl('^&', less_vec[j]))
              {
                less_vec[j] = gsub('(^[&])(.+)', '\\1!\\2', less_vec[j])
              } else
              {
                less_vec[j] = paste('!', less_vec[j], sep='')
              }
            }
            
          } else if(!grepl('!', more_vec[i]) & !grepl('!', less_vec[j])) #for both act terms. do not cross compare act and inh terms.
          {
            #Split & terms if present.
            mvec = unlist(strsplit(more_vec[i], '&'))
            lvec = unlist(strsplit(less_vec[j], '&'))
            
            #Check for presence of &.
            mand_bool = ifelse(grepl('&', more_vec[i]), T, F)
            land_bool = ifelse(grepl('&', less_vec[j]), T, F)
            
            mbool = mvec %in% lvec
            lbool = lvec %in% mvec
            
            #assign NA if there is no term remaining.
            if(length(mvec[!mbool])!=0)
            {
              if(!is.na(mvec[!mbool][1]))
              {
                more_vec[i] = paste(mvec[!mbool], collapse='&')
              } 
            } else
            {
              more_vec[i] = NA
            }
            
            if(length(lvec[!lbool])!=0)
            {
              if(!is.na(lvec[!lbool][1]))
              {
                less_vec[j] = paste(lvec[!lbool], collapse='&')
              } 
            } else
            {
              less_vec[j] = NA
            }
            
            #if this involves AND term, adds & indicator to term name.
            if(mand_bool & !is.na(more_vec[i]) & !grepl('&', more_vec[i]))
            {
              more_vec[i] = paste('&', more_vec[i], sep='')
            }
            if(land_bool & !is.na(less_vec[j]) & !grepl('&', less_vec[j]))
            {
              less_vec[j] = paste('&', less_vec[j], sep='')
            }
          }
        }
      }
    }
  }
  
  #Remove NAs.
  if(!is.null(more_vec))
  {
    more_vec = more_vec[!is.na(more_vec)]
  }
  if(!is.null(less_vec))
  {
    less_vec = less_vec[!is.na(less_vec)]
  }
  
  #If there are any AND term remaining here, they should be separated as they count as 2 steps.
  if(length(more_vec)!=0)
  {
    for(i in 1:length(more_vec))
    {
      if(grepl('.+&.+', more_vec[i]))
      {
        tmp_entry = unlist(strsplit(more_vec[i], '&'), use.names=F)
        names(tmp_entry) = rep(names(more_vec[i]), length(tmp_entry))
        more_vec = c(more_vec, tmp_entry)
        more_vec = more_vec[-i]
      }
    }
  }
  
  if(length(less_vec)!=0)
  {
    for(i in 1:length(less_vec))
    {
      if(grepl('.+&.+', less_vec[i]))
      {
        tmp_entry = unlist(strsplit(less_vec[i], '&'), use.names=F)
        names(tmp_entry) = rep(names(less_vec[i]), length(tmp_entry))
        less_vec = c(less_vec, tmp_entry)
        less_vec = less_vec[-i]
      }
    }
  }
  
  comb_vec = c(more_vec, less_vec)
  #Replace the gene names back into update functions.
  for(i in 1:length(x@target_var))
  {
    comb_vec = gsub(x@target_var[i], x@target[i], comb_vec)
  }
  
  stopifnot(all(!grepl('v[0-9]+s', unlist(strsplit(comb_vec, '&'))) | 
                  !grepl('^0$', unlist(strsplit(comb_vec, '&'))))) #check if all gene name conversions is successful.
  
  #Return results.
  if(directed)
  {
    return(list(more=more_vec, less=less_vec))
  } else
  {
    return(list(comb_vec))
  }
}