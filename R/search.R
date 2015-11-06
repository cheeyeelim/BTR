#' @title Training Model
#' 
#' @description
#' This function performs model training to find the best model, using information from data. It requires an initial state supplied to perform the search, and an initial model can also be supplied to be included in the initial population.
#' Note that if a model is supplied, and the genes in the model is different from the genes in the data, only the genes overlapping between model and data will be retained for further analysis.
#' 
#' @param bmodel Boolean model in data frame. If NULL, use a random Boolean model. Default to NULL.
#' @param edata list of 2 data frames. Initialised continuous and discretised expression data. Each data frame should have state(row) x gene(column).
#' @param istate data frame. Must have only 1 row, which represents 1 initial state.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' @param tol numeric. Specify the tolerance in ending condition. Default to 1e-6. It cannot be lower than .Machine$double.eps ^ 0.5.
#' @param inter_bool logical. Indicate whether to consider AND terms. Default to TRUE.
#' @param verbose logical. Specifies whether to give detailed output. Default to F.
#' @param self_loop logical. Indicates whether to allow self_loop in random starting model. Only used if is.null(bmodel). Default to F.
#' 
#' @export
model_train = function(bmodel=NULL, edata, istate, max_varperrule=6, tol=1e-6, inter_bool, verbose=F, self_loop=F)
{
  vcat('Preparing data for analysis.\n', verbose)
  
  if(class(edata)!='list' | length(edata)!=2)
  {
    stop('edata: Supply two expression data frames in a list.')
  }
  
  #Initialise input data.
  istate = initialise_data(istate, aslogic=T)
  cdata = initialise_data(edata[[1]])
  ddata = initialise_data(edata[[2]])
  
  #Initialise model.
  if(is.null(bmodel))
  {
    bmodel = gen_one_rmodel(colnames(istate), max_varperrule, inter_bool, self_loop)
  } else if(class(bmodel) != 'BoolModel')
  {
    bmodel = initialise_model(bmodel)
  }
  
  #Filtering expression data.
  stopifnot(colnames(edata[[1]])==colnames(edata[[2]]))
  overlap_gene = intersect(colnames(edata[[1]]), y=bmodel@target)
  nonoverlap_gene = bmodel@target[!(bmodel@target %in% overlap_gene)]
  names(overlap_gene) = bmodel@target_var[bmodel@target %in% overlap_gene]
  names(nonoverlap_gene) = bmodel@target_var[!(bmodel@target %in% overlap_gene)]

  fddata = filter_dflist(ddata, overlap_gene, F)
  fcdata = filter_dflist(cdata, overlap_gene, F)

  fcdata = unique_raw_data(fddata, fcdata) #removes duplicates in continuous data.
  fddata = unique(fddata)
  
  vcat('Start training.\n', verbose)

  #(3) Calling final combined search.
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
    
    #browser()
  }
  vcat(sprintf('Final iteration: %s.\n', cur_step), verbose)
  
  vcat('Stage 4: Performing consensus analysis.\n', verbose)
  consensus = model_consensus(best_model, inter_bool=inter_bool, max_varperrule=max_varperrule)
  
  output = list(consensus=consensus, best_model=best_model, best_score=best_score, ite_score=all_best_score, overlap_gene=overlap_gene, nonoverlap_gene=nonoverlap_gene)
  return(output)
}

#' @title Simplifying Model
#' 
#' @description
#' This method takes in a model and remove redundant terms wrt to a single initial state. 
#' Note that this model simplification is random, and the simplified model is not guaranteed to be the simplest model possible. It is only guaranteed to be a simpler model that can give the same state space as the orignal input model.
#' 
#' @param bmodel S4 BoolModel object.
#' @param istate data frame. Must have only 1 row, which represents 1 initial state.
#' @param inter_bool logical. Indicate whether to consider AND terms.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' @param verbose logical. Specifies whether to give detailed output. Default to F.
model_simplify = function(bmodel, istate, inter_bool, max_varperrule, verbose=F)
{
  vcat('Stage 1: Calculating score of initial model.\n', verbose)
  #Get the states of the original model.
  overlap_gene = bmodel@target
  fcdata = simulate_model(bmodel, istate)
  fcdata = fcdata+0 #convert logical to numeric.
  
  ori_score = calc_mscore(bmodel, istate, fcdata, overlap_gene, max_varperrule, simplify_bool=T)
  stopifnot(ori_score==0)
  
  next_bmodel = bmodel
  ite = 1
  while(TRUE)
  {
    cat(sprintf('Simplification iteration: %s\n', ite))
    #Generate list of minimally deleted models.
    
    vcat('Stage 2: Exploring neighbouring models.\n', verbose)
    
    mod_model = c(next_bmodel, minmod_model(next_bmodel, ibool=inter_bool, overlap_gene=overlap_gene)$del_list)
    #Breaking condition.
    if(length(mod_model) == 1) #model can no longer be simplified.
    {
      final_bmodel = mod_model[[1]]
      break
    }
    vcat(sprintf('Total neighbouring models: %s.\n', length(mod_model)), verbose)
    
    vcat('Stage 3: Simulating and calculating scores for models.\n', verbose)
    model_res = foreach(i=1:length(mod_model), .combine='c') %dopar% {
      model_score = calc_mscore(bmodel=mod_model[[i]], istate=istate, fcdata=fcdata, overlap_gene=overlap_gene, max_varperrule=max_varperrule, simplify_bool=T)
      names(model_score)=i
      model_score
    }
    
    all_final_score = unname(model_res)
    
    stopifnot(!is.null(model_res))
    stopifnot(length(all_final_score)==length(mod_model))
    
    #Breaking condition.
    if(!any(all_final_score[-1] == 0)) #model can no longer be simplified.
    {
      final_bmodel = mod_model[[1]]
      break
    }
    
    #Pick a random equivalent model for next iteration.
    best_ind = which.random.min(all_final_score)
    next_bmodel = mod_model[[best_ind]]
    
    ite = ite + 1
  }
  
  return(final_bmodel)
}

#' @title Intersection of input genes
#' 
#' @description
#' This function finds the intersection of input genes and provide a score for them. Return a consensus model or a vector of scores.
#' 
#' @param bmodel_list list of BoolModel.
#' @param inter_bool logical. Indicate whether to consider AND terms.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must not be smaller than number of variables. Default to 6.
#' @param format character. Specifies which format to return. Possible values: 'vec', 'df'. Default to 'vec'.
#' 
#' @export
model_consensus = function(bmodel_list, inter_bool, max_varperrule, format='vec')
{
  #(1) Convert all bmodels to encoded forms.
  encoding = get_encodings(bmodel_list[[1]], inter_bool)
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
