#' @title Get corresponding encodings for compression or decompression.
#' 
#' @description
#' This function gets all possible terms (single or double terms) and their corresponding encodings.
#' This function limits the number of possible variables in the model to 999.
#' 
#' @param bmodel S4 BoolModel object.
#' @param force_and logical. Whether to include ANDs in the encoding even if no AND is present in the bmodel object. Defaults to T.
#'
#' @export
get_encodings = function(bmodel, force_and=T)
{
  if(force_and)
  {
    and_bool = T
  } else
  {
    and_bool = check_and(bmodel)
  }
  
  
  #Get all possible terms.
  svar = bmodel@target_var
  if(and_bool)
  {
    dvar = sapply(combn(svar, 2, simplify=F), function(x) paste(x, collapse='&')) #get all possible interacting pairs.
    dvar = c(dvar, sapply(combn(svar, 2, simplify=F), function(x) paste(rev(x), collapse='&'))) #get the reversed pattern as well.
    
    term_pool = c(svar, dvar)
    term_pool = c('0', term_pool, '!0', paste('!', term_pool, sep='')) #add in inh terms.
    
    #Generate index for activation terms.
    num_pool = seq(1, length(svar)+1) #get numbers for svar.
    num_pool = c(num_pool, as.vector(replicate(2, seq(max(num_pool)+1, max(num_pool)+(length(dvar)/2))))) #get numbers for both forward and reverse dvar.
    
    #Generate index for inhibition terms.
    num_pool = c(num_pool, seq(max(num_pool)+1, max(num_pool)+length(svar)+1)) #get numbers for svar.
    num_pool = c(num_pool, as.vector(replicate(2, seq(max(num_pool)+1, max(num_pool)+(length(dvar)/2))))) #get numbers for both forward and reverse dvar.
    
    num_pool = unname(sapply(num_pool, function(x) ifelse(nchar(x)==1, paste('0', x, sep=''), x))) #convert single digit to double digit.
    num_pool = unname(sapply(num_pool, function(x) ifelse(nchar(x)==2, paste('0', x, sep=''), x))) #convert double digit to triple digit.
    num_pool = unname(sapply(num_pool, function(x) ifelse(nchar(x)==3, paste('0', x, sep=''), x))) #convert triple digit to quadruple digit.
    
    names(num_pool) = term_pool
  } else
  {
    term_pool = svar
    term_pool = c('0', term_pool, '!0', paste('!', term_pool, sep='')) #add in inh terms.
    
    #Generate index for activation terms.
    num_pool = seq(1, length(svar)+1) #get numbers for svar.
    
    #Generate index for inhibition terms.
    num_pool = c(num_pool, seq(max(num_pool)+1, max(num_pool)+length(svar)+1)) #get numbers for svar.
   
    num_pool = unname(sapply(num_pool, function(x) ifelse(nchar(x)==1, paste('0', x, sep=''), x))) #convert single digit to double digit.
    num_pool = unname(sapply(num_pool, function(x) ifelse(nchar(x)==2, paste('0', x, sep=''), x))) #convert double digit to triple digit.
    num_pool = unname(sapply(num_pool, function(x) ifelse(nchar(x)==3, paste('0', x, sep=''), x))) #convert triple digit to quadruple digit.
    
    names(num_pool) = term_pool
  }
  
  #Store original gene names.
  bnames = bmodel@target
  names(bnames) = bmodel@target_var

  stopifnot(all(!is.na(names(num_pool))))
  stopifnot(all(!is.na(term_pool)))
  
  return(list(num_pool, bnames))
}

#' @title Compress BoolModel
#' 
#' @description
#' This function compresses S4 BoolModel object by representing variables using numbers, and also only the act rules and inh rules are kept. 
#' Return a list of 3 vectors, corresponding to act rules and inh rules.
#' 
#' @param bmodel S4 BoolModel object.
#' @param encoding named numerical vector returned by get_encodings().
#'
#' @export
compress_bmodel = function(bmodel, encoding)
{
  encoding = encoding[[1]]
  stopifnot(length(bmodel@target)== length(encoding[!grepl('&|!', names(encoding))])-1) #check if bmodel and encoding match.
  
  #Convert variables into encodings.
  #Convert 'v1s' into '001'.
  act_list = lapply(bmodel@rule_act, function(x) unname(encoding[x]))
  inh_list = lapply(bmodel@rule_inh, function(x) unname(encoding[paste('!', x, sep='')]))
  
  stopifnot(all(!is.na(unlist(act_list))))
  stopifnot(all(!is.na(unlist(inh_list))))
  
  #Get max number of variables.
  max_lim = max(sapply(act_list, length) + sapply(inh_list, length))
  
  #Convert '0001' into 1.0001.
  act_list = lapply(1:length(act_list), function(x) as.numeric(paste(x, act_list[[x]], sep='.')))
  inh_list = lapply(1:length(inh_list), function(x) as.numeric(paste(x, inh_list[[x]], sep='.')))
  
  #Convert list into vector.
  act_vec = unique(unlist(act_list))
  inh_vec = unique(unlist(inh_list))
  
  #Add in integer that represents empty spots, according to max_lim.
  empty_term = max_lim - (rle(round(act_vec))$lengths+rle(round(inh_vec))$lengths) #get freq of missing terms.
  empty_term[empty_term < 0] = 0 #set negative values to 0. note that negative values can occur because this is a soft constraint.
  empty_term = unlist(sapply(1:length(empty_term), function(x) rep(x, empty_term[x])))
  
  return(sort(c(act_vec,inh_vec, empty_term)))
}

#' @title Decompress BoolModel
#' 
#' @description
#' This function decompresses the bmodel compressed by compress_bmodel(). 
#' Return a S4 BoolModel object.
#' 
#' @param x vector returned by compress_bmodel.
#' @param encoding named numerical vector returned by get_encodings().
#' @param format character. Specifies which format to return. Possible values: 'bmodel', 'df', 'amat', 'simp_df'. Default to 'bmodel'.
#'
#' @export
decompress_bmodel = function(x, encoding, format='bmodel')
{
  gene = encoding[[2]]
  encoding = encoding[[1]]
  
  #Setup default gene names.
  gene_var = names(encoding)[!(grepl('&',names(encoding))|grepl('!',names(encoding)))][-1]
  if(is.null(gene))
  {
    gene = gene_var
  }
  stopifnot(length(gene_var)==length(gene))
  
  #Round up values to 4 decimals. (due to external solver returning continuous values.)
  x = round(x, 4)
  
  #Remove all empty terms from x.
  x = x[!(x%%1==0)] #check for integer.
  x = x[!(sapply(x, function(y) isTRUE(all.equal(y%%1,0.9999))))] #check for lower_bound values (i.e. 2.9999). (it is equivalent to integer.)
  
  #Filter encodings to remove reversed AND terms, if present.
  encoding = encoding[!duplicated(encoding)]
  
  #Convert encodings into variables.
  #Obtain 2 from 2.001.
  comb_ind = round(x)
  
  #Convert numbers to characters.
  comb_rule = as.character(x)
  
  #Convert '2' to '2.'.
  comb_rule = unname(sapply(comb_rule, function(y) ifelse(!grepl('[.]', y), paste(y, '.', sep=''), y)))
  
  #Convert '2.' to '2.0'.
  comb_rule = unname(sapply(comb_rule, function(y) ifelse(grepl('[.]$', y), paste(y, '0', sep=''), y)))
  
  #Convert '2.0' to '2.00'.
  comb_rule = unname(sapply(comb_rule, function(y) ifelse(grepl('[.].$', y), paste(y, '0', sep=''), y)))
 
  #Convert '2.00' to '2.000'.
  comb_rule = unname(sapply(comb_rule, function(y) ifelse(grepl('[.]..$', y), paste(y, '0', sep=''), y)))
  
  #Convert '2.000' to '2.0000'.
  comb_rule = unname(sapply(comb_rule, function(y) ifelse(grepl('[.]...$', y), paste(y, '0', sep=''), y)))
  
  #Obtain '0001' from '2.0001'.
  comb_rule = gsub('^[0-9]+[.]([0-9]+)', '\\1', comb_rule)
  
  #print(comb_rule)
  stopifnot(all(nchar(comb_rule)==4))
  
  #Remove 0 and !0.
  zero_enc = encoding[grepl('0$', names(encoding))] #get encodings for 0 and !0.
  stopifnot(length(zero_enc)==2)
  comb_ind = comb_ind[comb_rule!=zero_enc[1]] #remove 0.
  comb_rule = comb_rule[comb_rule!=zero_enc[1]] #remove 0.
  comb_ind = comb_ind[comb_rule!=zero_enc[2]] #remove !0.
  comb_rule = comb_rule[comb_rule!=zero_enc[2]] #remove !0.
  stopifnot(length(comb_ind)==length(comb_rule))
  
  #Generate reverse map.
  rev_encoding = names(encoding)
  names(rev_encoding) = encoding
  
  #Main conversion.
  comb_rule = unname(rev_encoding[comb_rule])
  stopifnot(all(!is.na(comb_rule)))
  
  if(format=='df')
  {
    #Converting rules from var to gene names.
    conv_comb_rule = comb_rule
    for(i in 1:length(gene_var))
    {
      conv_comb_rule = gsub(gene_var[i], gene[i], conv_comb_rule)
    }
    conv_comb_rule = unname(split(conv_comb_rule, as.factor(comb_ind)))
    
    #If there is completely empty rules for a gene, add placeholder in.
    if(length(unique(comb_ind)) != length(gene_var))
    {
      missing_ind = which(!(1:length(gene_var) %in% unique(comb_ind)))
      
      for(i in missing_ind)
      {
        conv_comb_rule = append(conv_comb_rule, '0', i-1)
      }
    }
    stopifnot(length(conv_comb_rule)==length(gene_var))
  }
  
  #Convert vector into list.
  comb_rule = unname(split(comb_rule, as.factor(comb_ind)))

  #If there is completely empty rules for a gene, add placeholder in.
  if(length(unique(comb_ind)) != length(gene_var))
  {
    missing_ind = which(!(1:length(gene_var) %in% unique(comb_ind)))
    
    for(i in missing_ind)
    {
      comb_rule = append(comb_rule, '0', i-1)
    }
  }
  stopifnot(length(comb_rule)==length(gene_var))
  
  if(format=='bmodel')
  {
    #Split into act and inh list.
    act_rule = lapply(comb_rule, function(x) x[!grepl('!', x)])
    inh_rule = lapply(comb_rule, function(x) x[grepl('!', x)])
    
    #Set '0' for empty rules.
    act_rule[which(sapply(act_rule, function(x) length(x))==0)] = list('0')
    inh_rule[which(sapply(inh_rule, function(x) length(x))==0)] = list('0')
    
    inh_rule = lapply(inh_rule, function(x) gsub('!', '', x)) #remove !.
    
    stopifnot(length(gene)==length(act_rule))
    
    output = BoolModel(target=gene, target_var=gene_var, 
                           rule_act=act_rule, rule_inh=inh_rule)
  } else if(format=='df')
  {
    #This way of preparing output ensures that the order of 'output' is the same as 'x'.
    out_var = c()
    out_gene = c()
    for(i in 1:length(comb_rule))
    {
      for(j in 1:length(comb_rule[[i]]))
      {
        if(!grepl('!', comb_rule[[i]][j]))
        {
          if(comb_rule[[i]][j] != '0')
          {
            out_var = c(out_var, paste(gene_var[i], '<-', comb_rule[[i]][j]))
            out_gene = c(out_gene, paste(gene[i], '<-', conv_comb_rule[[i]][j]))
          }
        } else
        {
          if(comb_rule[[i]][j] != '0')
          {
            out_var = c(out_var, paste(gene_var[i], '|-', gsub('!', '', comb_rule[[i]][j])))
            out_gene = c(out_gene, paste(gene[i], '|-', gsub('!', '', conv_comb_rule[[i]][j])))
          }
        }
      }
    }
    
    output = cbind(out_gene, out_var)
  } else if(format=='amat')
  {
    comb_rule = lapply(1:length(comb_rule), function(x) unlist(strsplit(comb_rule[[x]], '&'))) #remove &.
    comb_rule = lapply(1:length(comb_rule), function(x) unique(gsub('!', '', comb_rule[[x]]))) #remove ! and non-unique gene interactions.
    comb_rule = lapply(1:length(comb_rule), function(x) comb_rule[[x]][order(as.numeric(gsub('v(.+)s', '\\1', comb_rule[[x]])))]) #reorder the elements.
    
    output = matrix(NA, ncol=length(gene_var), nrow=length(gene_var))
    for(i in 1:nrow(output))
    {
      output[i,] = gene_var %in% comb_rule[[i]]
    }
    output = t(output)
    output = output + 0
    
    colnames(output) = gene
    rownames(output) = gene
  } else if(format=='direct_amat')
  {
    comb_rule = lapply(1:length(comb_rule), function(x) unlist(strsplit(comb_rule[[x]], '&'))) #remove &.
    comb_rule = lapply(comb_rule, unique) #remove non-unique gene interactions.
    
    output = matrix(0, ncol=length(gene_var), nrow=length(gene_var))
    for(i in 1:length(comb_rule))
    {
      for(j in 1:length(comb_rule[[i]]))
      {
        row_ind = as.numeric(gsub('.+([0-9]+)s', '\\1', comb_rule[[i]][j]))
        if(grepl('!', comb_rule[[i]][j]))
        {
          output[row_ind,i] = -1
        } else
        {
          output[row_ind,i] = 1
        }
      }
    }
    
    colnames(output) = gene
    rownames(output) = gene
  } else if(format=='simp_df')
  {
    comb_rule = lapply(1:length(comb_rule), function(x) unlist(strsplit(comb_rule[[x]], '&'))) #remove &.
    comb_rule = lapply(1:length(comb_rule), function(x) unique(gsub('!', '', comb_rule[[x]]))) #remove ! and non-unique gene interactions.
    comb_rule = lapply(1:length(comb_rule), function(x) comb_rule[[x]][order(as.numeric(gsub('v(.+)s', '\\1', comb_rule[[x]])))]) #reorder the elements.
    comb_rule = comb_rule[!sapply(comb_rule, function(x) x[1]=='0')] #remove empty rules.
    
    output = matrix(NA, ncol=2, nrow=length(unlist(comb_rule)))
    colnames(output) = c('from', 'to')
    
    entry_ind = 1
    for(i in 1:length(comb_rule))
    {
      for(j in 1:length(comb_rule[[i]]))
      {
        output[entry_ind,1] = comb_rule[[i]][j]
        output[entry_ind,2] = paste('v', i, 's', sep='')
        
        entry_ind = entry_ind + 1
      }
    }
  } else
  {
    stop('Please provide a valid format.')
  }
  
  return(output)
}
