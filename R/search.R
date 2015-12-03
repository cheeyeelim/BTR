#' @title Training Model
#' 
#' @description
#' This function performs model training to find the best model, using information from data. It requires an initial state supplied to perform the search, and an initial model can also be supplied to be included in the initial population.
#' Note that if a model is supplied, and the genes in the model is different from the genes in the data, only the genes overlapping between model and data will be retained for further analysis.
#' 
#' @param cdata data frame of expression data. Should have state(row) x gene(column).
#' @param ddata discretised data frame of expression data. Must supply when preprocess=F. Obtain from initialise_raw_data(). Defaults to NULL.
#' @param bmodel Boolean model in data frame. If NULL, use a random Boolean model. Defaults to NULL.
#' @param istate data frame. Must have only 1 row, which represents 1 initial state. Defaults to NULL.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must be higher than number of genes. Defaults to 6.
#' @param and_bool logical. Whether to consider AND terms. IF bmodel is not NULL, defaults to whether AND interaction is included in bmodel. If bmodel is NULL, then defaults to TRUE.
#' @param self_loop logical. Whether to allow self_loop in random starting model. Only used if is.null(bmodel). Default to F.
#' @param con_thre numerical. Threshold used to generating the final consensus model. Must be between 0 and 1.
#' @param tol numeric. Tolerance in ending condition. Default to 1e-6. It cannot be lower than .Machine$double.eps ^ 0.5.
#' @param verbose logical. Whether to give detailed output to the screen. Defaults to F.
#' @param detailed_output logical. Whether to return only the model inferred, or all the details obtained during optimisation. Defaults to F.
#' 
#' @export
model_train = function(cdata, ddata=NULL, bmodel=NULL, istate=NULL, max_varperrule=6, and_bool=T, self_loop=F, con_thre=0.3, tol=1e-6, verbose=F, detailed_output=F)
{
  vcat('Preparing data for analysis.\n', verbose)
  
  #Initialise model.
  if(is.null(bmodel))
  {
    bmodel = gen_one_rmodel(colnames(cdata), max_varperrule, and_bool, self_loop)
  } else
  {
    if(class(bmodel) != 'BoolModel')
    {
      bmodel = initialise_model(bmodel)
    }
    
    if(check_and(bmodel) != and_bool)
    {
      and_bool = check_and(bmodel)
    }
  }
  
  #Initialise initial state.
  if(is.null(istate))
  {
    istate = rbinom(length(bmodel@target), 1, 0.5)
    #Getting a random initial state.
    while(mean(istate) > 0.9 | mean(istate) < 0.1) #do not want initial state that is too homogenous.
    {
      istate = rbinom(length(bmodel@target), 1, 0.5)
    }
    istate = data.frame(matrix(istate, nrow=1))
    colnames(istate) = bmodel@target
  }
  istate = initialise_data(istate, aslogic=T)
  
  #Filtering expression data.
  overlap_gene = intersect(colnames(cdata), y=bmodel@target)
  nonoverlap_gene = bmodel@target[!(bmodel@target %in% overlap_gene)]
  names(overlap_gene) = bmodel@target_var[bmodel@target %in% overlap_gene]
  names(nonoverlap_gene) = bmodel@target_var[!(bmodel@target %in% overlap_gene)]
  
  fddata = filter_dflist(ddata, overlap_gene, F)
  fcdata = filter_dflist(cdata, overlap_gene, F)
  
  fcdata = unique_raw_data(fddata, fcdata) #removes duplicates in continuous data.
  fddata = unique(fddata)
  
  vcat('Start training.\n', verbose)
  
  #(3) Calling final combined search.
  i = 0 #suppress check error on non-visible global binding.
  best_model = c()
  best_score = c()
  all_best_score = list()
  cur_step = 1
  next_break_ite = 1
  while(TRUE)
  {
    vcat(sprintf('Current iteration: %s.\n', cur_step), verbose)
    
    if(cur_step == 1)
    {
      if(length(bmodel)==1)
      {
        mod_model = list(bmodel)
      } else
      {
        mod_model = bmodel
      }
    } else
    {
      mod_model = next_model
    }
    
    vcat('Stage 1: Exploring neighbouring models.\n', verbose)
    mod_model = foreach(i=1:length(mod_model)) %dopar% {
      c(mod_model[[i]], unlist(minmod_model(mod_model[[i]], overlap_gene=overlap_gene)))
    }
    mod_model = unlist(mod_model)
    
    vcat(sprintf('Total neighbouring models: %s.\n', length(mod_model)), verbose)
    if(length(mod_model)>1000000)
    {
      #due to memory limit and runtime consideration.
      break
    }
    
    vcat('Stage 2: Simulating and calculating scores for models.\n', verbose)
    all_final_score = foreach(i=1:length(mod_model)) %dopar% {
      model_score = calc_mscore(bmodel=mod_model[[i]], istate=istate, fcdata=fcdata, overlap_gene=overlap_gene, max_varperrule=max_varperrule)
      unname(model_score['f'])
    }
    
    stopifnot(all(sapply(all_final_score, function(x) !is.null(x)))) #there will be NULL if any problem occurs within the workers of foreach.
    stopifnot(all(sapply(all_final_score, length) == sapply(mod_model, length)))
    
    all_final_score = unlist(all_final_score)
    
    #Locating all equally best score branches.
    best_ind = which(all_final_score == min(all_final_score)) #take all min scores.
    best_score = all_final_score[best_ind]
    best_model = mod_model[best_ind]
    
    next_model = best_model
    
    vcat(sprintf('Total model retained for next iteration: %s.\n', length(next_model)), verbose)
    
    if(cur_step > 1) #checking convergence, and breaking condition.
    {
      if(any(best_score == 0)) #if meet any model with 0 score,  #if the score changes less than the tolerance between 2 iterations.
      {
        break
      } else if(isTRUE(all.equal(mean(previous_score),mean(best_score), tolerance=tol)))
      {
        #increase break count. Iteration will be broken once the count reaches 2.
        next_break_ite = next_break_ite + 1
        if(next_break_ite == 2)
        {
          break
        }
      }  else if(length(next_model) == 0)
      {
        stop('Error: length(next_model)==0.')
      } else
      {
        next_break_ite = 1
      }
    } 
    
    previous_score = best_score #store it for comparison.
    all_best_score = c(all_best_score, list(best_score))
    cur_step = cur_step + 1
  }
  vcat(sprintf('Final iteration: %s.\n', cur_step), verbose)
  
  vcat('Stage 4: Performing consensus analysis.\n', verbose)
  consensus = model_consensus(best_model, max_varperrule=max_varperrule)
  res_con = consensus[consensus > con_thre]
  final_model = decompress_bmodel(as.numeric(names(res_con)), get_encodings(bmodel), gene=bmodel@target)
  
  if(detailed_output)
  {
    output = list(consensus=consensus, final_model=final_model, best_model=best_model, 
                  best_score=best_score, ite_score=all_best_score, overlap_gene=overlap_gene, nonoverlap_gene=nonoverlap_gene)
  } else
  {
    output = final_model
  }
  
  return(output)
}

#' @title Intersection of genes
#' 
#' @description
#' This function finds the intersection of genes and provide a score for them. Return a consensus model or a vector of scores.
#' 
#' @param bmodel_list list of BoolModel.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' @param format character. Specifies which format to return. Possible values: 'vec', 'df'. Default to 'vec'.
model_consensus = function(bmodel_list, max_varperrule, format='vec')
{
  #(1) Convert all bmodels to encoded forms.
  encoding = get_encodings(bmodel_list[[1]])
  encmodel_list = lapply(bmodel_list, function(x) compress_bmodel(x, encoding, max_varperrule))
  
  #(2) Check and remove duplicated models.
  #stopifnot(all(sapply(encmodel_list, length)==length(encmodel_list[[1]]))) #all models must have same lengths after encoding.
  encmodel_list = unique(encmodel_list)
  
  #(3) Count frequency of each term.
  rle_freq = rle(sort(unlist(encmodel_list)))
  term_freq = rle_freq$lengths
  names(term_freq) = rle_freq$values
  
  #(4) Remove integers (i.e. empty terms)
  term_freq = term_freq[as.numeric(names(term_freq))%%1!=0]
  
  #(5) Remove 0 and !0.
  zero_enc = encoding[grepl('0$', names(encoding))] #get encodings for 0 and !0.
  stopifnot(length(zero_enc)==2)
  term_freq = term_freq[!grepl(paste(zero_enc[1], '$', sep=''), names(term_freq))] #remove 0.
  term_freq = term_freq[!grepl(paste(zero_enc[2], '$', sep=''), names(term_freq))] #remove !0.
  
  #(6) Obtain final frequency.
  final_res = signif(term_freq/length(encmodel_list)) #convert into percentage.
  
  if(format=='df')
  {
    #(7) Decompress the genes into a Boolean model.
    out_bmodel = decompress_bmodel(as.numeric(names(term_freq)), encoding, gene=bmodel_list[[1]]@target, format='df')
    stopifnot(length(term_freq)==nrow(out_bmodel))
    
    output = cbind(out_bmodel, final_res)
  } else if(format=='vec')
  {
    output = final_res
  } else
  {
    stop('Invalid format parameter.')
  }

  return(output)
}
