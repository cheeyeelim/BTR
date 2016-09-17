#' @title Training Model
#' 
#' @description
#' This function performs model training to find the best model, using information from data. It requires an initial state supplied to perform the search, and an initial model can also be supplied to be included in the initial population.
#' Note that if a model is supplied, and the genes in the model is different from the genes in the data, only the genes overlapping between model and data will be retained for further analysis.
#' 
#' @param cdata data frame of expression data. Should have state(row) x gene(column).
#' @param bmodel Boolean model in data frame. If NULL, use a random Boolean model. Defaults to NULL.
#' @param istate data frame. Must have only 1 row, which represents 1 initial state. Defaults to NULL.
#' @param max_varperrule integer. Maximum number of terms per rule (combining both act and inh rule). Note that this number must be higher than number of genes. Defaults to 3.
#' @param and_bool logical. Whether to consider AND terms. IF bmodel is not NULL, defaults to whether AND interaction is included in bmodel. If bmodel is NULL, then defaults to TRUE.
#' @param self_loop logical. Whether to allow self_loop in random starting model. Default to F.
#' @param con_thre numerical. Threshold used to generating the final consensus model. Must be between 0 and 1.
#' @param tol numeric. Tolerance in ending condition. Default to 1e-6. It cannot be lower than .Machine$double.eps ^ 0.5.
#' @param verbose logical. Whether to give detailed output to the screen. Defaults to F.
#' @param detailed_output logical. Whether to return only the model inferred, or all the details obtained during optimisation. Defaults to F.
#' 
#' @examples
#' data(wilson_raw_data)
#' cdata = initialise_raw_data(wilson_raw_data, max_expr = 'low')
#' 
#' #select only relevant cells.
#' cell_ind = grepl('cmp', rownames(cdata)) | grepl('gmp', rownames(cdata)) 
#' fcdata = cdata[cell_ind,]
#' 
#' #select genes to be included.
#' gene_ind = c('fli1', 'gata1', 'gata2', 'gfi1', 'scl', 'sfpi1') 
#' fcdata = fcdata[, gene_ind]
#' 
#' final_model = model_train(cdata=fcdata, max_varperrule=2)
#' plotBM(final_model)
#' 
#' @export
model_train = function(cdata, bmodel=NULL, istate=NULL, max_varperrule=6, and_bool=T, self_loop=F, con_thre=0.3, tol=1e-6, verbose=F, detailed_output=F)
{
  vcat('Preparing data for analysis.\n', verbose)
  
  #Initialise model.
  if(is.null(bmodel))
  {
    bmodel = gen_one_rmodel(var=colnames(cdata), mvar=max_varperrule, and_bool=and_bool, self_loop=self_loop)
  } else
  {
    if(class(bmodel) != 'BoolModel')
    {
      bmodel = initialise_model(bmodel)
    }
  }
  
  if(length(bmodel)==1)
  {
    tmp_bm = bmodel
  } else
  {
    tmp_bm = bmodel[[1]]
  }
  
  #Initialise initial state.
  if(is.null(istate))
  {
    istate = rbinom(length(tmp_bm@target), 1, 0.5)
    #Getting a random initial state.
    while(mean(istate) > 0.9 | mean(istate) < 0.1) #do not want initial state that is too homogenous.
    {
      istate = rbinom(length(tmp_bm@target), 1, 0.5)
    }
    istate = data.frame(matrix(istate, nrow=1))
    colnames(istate) = tmp_bm@target
  }
  istate = initialise_data(istate, aslogic=T)
  
  #Filtering expression data.
  overlap_gene = intersect(colnames(cdata), y=tmp_bm@target)
  nonoverlap_gene = tmp_bm@target[!(tmp_bm@target %in% overlap_gene)]
  names(overlap_gene) = tmp_bm@target_var[tmp_bm@target %in% overlap_gene]
  names(nonoverlap_gene) = tmp_bm@target_var[!(tmp_bm@target %in% overlap_gene)]
  
  fcdata = filter_dflist(cdata, overlap_gene, F)
  
  model_encoding = get_encodings(tmp_bm) #for compression purpose.
  
  vcat('Start training.\n', verbose)
  
  #(3) Calling final combined search.
  i = 0 #suppress check error on non-visible global binding.
  #all_best_score = list()
  start_step = 1
  max_step = 100
  best_model = c()
  best_score = c()
  next_break_ite = 1
  for(cur_step in start_step:max_step)
  {
    vcat(sprintf('Current iteration: %s.\n', cur_step), verbose)

    if(cur_step == 1)
    {
      if(length(bmodel)==1)
      {
        mod_model = list(compress_bmodel(bmodel, model_encoding))
      } else
      {
        mod_model = lapply(bmodel, function(x) compress_bmodel(x, model_encoding))
      }
    } else
    {
      mod_model = next_model
    }
    
    vcat('Stage 1: Exploring neighbouring models.\n', verbose)
    mod_model = foreach(i=1:length(mod_model), .combine = c) %dopar% {
      c(list(mod_model[[i]]), minmod_model(mod_model[[i]], overlap_gene=overlap_gene, encoded = T, model_encoding=model_encoding, and_bool=and_bool, self_loop=self_loop))
    }
    stopifnot(!any(is.null(mod_model))) #there will be NULL if any problem occurs within the workers of foreach.
    
    mod_model = unique(mod_model)
    
    vcat(sprintf('Total neighbouring models: %s.\n', length(mod_model)), verbose)

    #Take random 1e6 models, if total number of neighbours more than that.
    if(length(mod_model)>1000000)
    {
      mod_model = sample(mod_model, 1000000)
    }
    
    vcat('Stage 2: Simulating and calculating scores for models.\n', verbose)
    all_final_score = foreach(i=1:length(mod_model), .inorder=F, .combine = c) %dopar% {
      model_score = calc_mscore(bmodel=mod_model[[i]], istate=istate, fcdata=fcdata, overlap_gene=overlap_gene, max_varperrule=max_varperrule, model_encoding=model_encoding)
      unname(model_score['f'])
    }

    stopifnot(!any(is.null(all_final_score))) #there will be NULL if any problem occurs within the workers of foreach.
    stopifnot(length(all_final_score) == length(mod_model))
    
    #Locating all equally best score branches.
    best_ind = which(all_final_score == min(all_final_score)) #take all min scores.
    best_score = all_final_score[best_ind]
    best_model = mod_model[best_ind]
    #save(best_score, best_model, file='tmp_cur_best_result.rda')
    
    next_model = unique(best_model)
    vcat(sprintf('Current mean best score: %s.\n', round(mean(best_score), 6)), verbose)
    
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
    #all_best_score = c(all_best_score, list(best_score))

    gc()
  }
  vcat(sprintf('Final iteration: %s.\n', cur_step), verbose)
  
  vcat('Stage 4: Performing consensus analysis.\n', verbose)
  consensus = model_consensus(best_model, model_encoding=model_encoding)
  res_con = consensus[consensus > con_thre]
  final_model = decompress_bmodel(as.numeric(names(res_con)), model_encoding)
  
  if(detailed_output)
  {
      output = list(consensus=consensus, final_model=final_model, best_model=best_model, 
                    best_score=best_score, overlap_gene=overlap_gene, nonoverlap_gene=nonoverlap_gene)
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
#' @param model_encoding list. Use for compressing and decompressing Boolean models.
#' @param format character. Specifies which format to return. Possible values: 'vec', 'df'. Default to 'vec'.
model_consensus = function(bmodel_list, model_encoding, format='vec')
{
  model_encoding = model_encoding[[1]]
  if(class(bmodel_list[[1]])=='BoolModel')
  {
    #(1) Convert all bmodels to encoded forms.
    encmodel_list = lapply(bmodel_list, function(x) compress_bmodel(x, model_encoding))
  } else
  {
    encmodel_list = bmodel_list
  }
  
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
  zero_enc = model_encoding[grepl('0$', names(model_encoding))] #get model_encodings for 0 and !0.
  stopifnot(length(zero_enc)==2)
  term_freq = term_freq[!grepl(paste(zero_enc[1], '$', sep=''), names(term_freq))] #remove 0.
  term_freq = term_freq[!grepl(paste(zero_enc[2], '$', sep=''), names(term_freq))] #remove !0.
  
  #(6) Obtain final frequency.
  final_res = signif(term_freq/length(encmodel_list)) #convert into percentage.
  
  if(format=='df')
  {
    #(7) Decompress the genes into a Boolean model.
    out_bmodel = decompress_bmodel(as.numeric(names(term_freq)), model_encoding, format='df')
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
