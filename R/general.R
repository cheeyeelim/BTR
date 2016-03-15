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

#' @title Check if containing AND terms
#' 
#' @description
#' This function checks if a particular Boolean model contains AND terms.
#' 
#' @param bmodel BoolModel object.
check_and = function(bmodel)
{
  and_bool = F
  
  if(any(grepl('&', bmodel@rule_act)))
  {
    and_bool = T
  } else if(any(grepl('&', bmodel@rule_inh)))
  {
    and_bool = T
  }
  
  return(and_bool)
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

# #' @title Filter a Boolean model
# #' 
# #' @description
# #' This function filters a Boolean model to retain only overlapping genes.
# #' 
# #' @param x BoolModel.
# #' @param y character vector. Overlapping genes.
# filter_bmodel = function(x, y)
# {
#   stopifnot(class(x)== 'BoolModel')
#   
#   if(uniq_bool)
#   {
#     return(unique(x[,colnames(x) %in% y, drop=F]))
#   } else
#   {
#     return(x[,colnames(x) %in% y, drop=F])
#   }
# }

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

#' @title Calculate distance between Boolean models
#' 
#' @description
#' This method takes in two models and calculate the distance between them. The value return indicate the number of steps between the two models.
#' 
#' @param x S4 BoolModel object. Test model.
#' @param y S4 BoolModel object. Reference model.
#' 
#' @export
model_dist = function(x, y)
{
  set_diff = unlist(model_setdiff(x, y))
  
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
#' @param directed logical. If TRUE, return the difference in terms with respect to x.
#' 
#' @export
model_setdiff = function(x, y, directed=F)
{
  stopifnot(length(x@target) == length(y@target))
  stopifnot(get_encodings(x)[[1]]==get_encodings(y)[[1]])
  
  ind = get_encodings(x)
  
  x1 = compress_bmodel(x, ind)
  x2 = compress_bmodel(y, ind)
  
  ind = ind[[1]]
  
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