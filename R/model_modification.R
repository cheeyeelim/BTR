#' @title Minimal modification of whole Boolean model
#' 
#' @description
#' This function generates all possible boolean models minimally modified. Returns a lists of 2 lists, deleted models and added models
#' 
#' @param bm S4 BoolModel object.
#' @param index integer. Specifying rule of which gene to modify. If NULL, modifies all rules in the model. Defaults to NULL.
#' @param overlap_gene character vector. Specify which genes are present in both model and data inputs.
#' @param sep_list logical. Separate add and del lists.
#' @param encoded logical. Return Boolean models in encoded form to save space.
#' @param model_encoding list. Only used if encoded=T.
#' @param and_bool logical. Default to check_and().
#' @param self_loop logical. Whether to allow self_loop in random starting model. Default to F.
#' 
#' @export
minmod_model = function(bm, index=NULL, overlap_gene=NULL, sep_list=F, encoded=F, model_encoding, and_bool=NULL, self_loop=F)
{
  if(class(bm) != 'BoolModel')
  {
    bm = decompress_bmodel(bm, model_encoding)
  }
  
  if(is.null(index))
  {
    #Modify only the overlapping genes
    gene_ind = which(bm@target %in% overlap_gene)
    
    del_list = c()
    add_list = c()
    for(ind in gene_ind) #Modify act and inh rules for one gene at a time.
    {
      tmp_list = minmod_internal(bm, ind, encoded, model_encoding, and_bool, self_loop)
      if(!is.null(tmp_list$dellist))
      {
        del_list = c(del_list, tmp_list$dellist)
      }
      if(!is.null(tmp_list$addlist))
      {
        add_list = c(add_list, tmp_list$addlist)
      }
    }
    
    if(sep_list)
    {
      out_list = list(del_list=del_list,
                      add_list=add_list)
    } else
    {
      out_list = c(del_list, add_list)
    }
  } else
  {
    if(sep_list)
    {
      out_list = minmod_internal(bm, index, encoded, model_encoding, and_bool, self_loop)
    } else
    {
      out_list = unlist(minmod_internal(bm, index, encoded, model_encoding, and_bool, self_loop))
    }
  }

  return(out_list)
}

#' @title Inner function of minimal modification of whole Boolean model
#' 
#' @description
#' This function generates all possible boolean models minimally modified. Returns a lists of 2 lists, deleted models and added models
#' 
#' @param bm S4 BoolModel object.
#' @param index integer. Specifying rule of which gene to modify.
#' @param encoded logical. Return Boolean models in encoded form to save space.
#' @param model_encoding list. Only used if encoded=T.
#' @param and_bool logical. Default to check_and() if unspecified.
#' @param self_loop logical. Whether to allow self_loop in random starting model. Default to F.
minmod_internal = function(bm, index, encoded, model_encoding, and_bool, self_loop=F)
{
  if(is.null(and_bool))
  {
    and_bool = check_and(bm)
  }
  
  arule = bm@rule_act[[index]]
  irule = bm@rule_inh[[index]]
  
  #Deletion of act rule.
  dellist_arule = list()
  if(length(arule) == 1)
  {
    if(grepl('&' , arule))
    {
      tmp = unlist(strsplit(arule[1], '&')) #split the AND term into 2.
      dellist_arule = c(dellist_arule, list(c(tmp[1], arule[-1]), c(tmp[2], arule[-1])))
    } else  #allow empty rule.
    {
      dellist_arule = c(dellist_arule, list(arule[-1]))
    }
#     else if(arule[1] != '0' & irule[1] != '0')  #not allowing empty rule.
#     {
#       dellist_arule = c(dellist_arule, list('0'))
#     }
  } else
  {
    for(i in 1:length(arule))
    {
      if(grepl('&', arule[i]))
      {
        tmp = unlist(strsplit(arule[i], '&')) #split the AND term into 2.
        dellist_arule = c(dellist_arule, list(c(tmp[1], arule[-i]), c(tmp[2], arule[-i])))
      } else
      {
        dellist_arule = c(dellist_arule, list(arule[-i]))
      }
    }
  }
  dellist_arule[sapply(dellist_arule, length)==0] = list('0') #if there is any empty term at the end, add in '0'

  #Deletion of inh rule.
  dellist_irule = list()
  if(length(irule) == 1)
  {
    if(grepl('&' , irule))
    {
      tmp = unlist(strsplit(irule[1], '&')) #split the AND term into 2.
      dellist_irule = c(dellist_irule, list(c(tmp[1], irule[-1]), c(tmp[2], irule[-1])))
    } else  #allow empty rule.
    {
      dellist_irule = c(dellist_irule, list(irule[-1]))
    }
#     } else if(arule[1] != '0' & irule[1] != '0') #not allowing empty rule.
#     {
#       dellist_irule = c(dellist_irule, list('0'))
#     }
  } else
  {
    for(i in 1:length(irule))
    {
      if(grepl('&', irule[i]))
      {
        tmp = unlist(strsplit(irule[i], '&')) #split the AND term into 2.
        dellist_irule = c(dellist_irule, list(c(tmp[1], irule[-i]), c(tmp[2], irule[-i])))
      } else
      {
        dellist_irule = c(dellist_irule, list(irule[-i]))
      }
    }
  }
  dellist_irule[sapply(dellist_irule, length)==0] = list('0') #if there is any empty term at the end, add in '0'
  
  #Addition of act rule. (single)
  if(self_loop)
  {
    pos_actterm = bm@target_var[!(bm@target_var %in% unlist(strsplit(c(arule, irule), '&')))]
  } else
  {
    pos_actterm = bm@target_var[!(bm@target_var %in% unlist(strsplit(c(arule, irule), '&')) | bm@target_var %in% bm@target_var[index])]
  }

  addlist_arule = list()
  if(length(pos_actterm) != 0)
  {
    if(arule[1] == '0')
    {
      for(i in 1:length(pos_actterm))
      {
        addlist_arule = c(addlist_arule, list(pos_actterm[i]))
      }
    } else
    {
      for(i in 1:length(pos_actterm))
      {
        addlist_arule = c(addlist_arule, list(c(arule, pos_actterm[i])))
      }
    }
    
    if(and_bool)
    {
      #Addition of act rule. (double)
      if(arule[1] != '0')
      {
        for(i in 1:length(pos_actterm))
        {
          for(j in 1:length(arule))
          {
            if(grepl('&', arule[j]))
            {
              next
            } else
            {
              tmp = sprintf('%s&%s', pos_actterm[i], arule[j])
              addlist_arule = c(addlist_arule, list(c(arule[-j], tmp)))
            }
          }
        }
      }
    }
  }
  
  #Addition of inh rule. (single)
  if(self_loop)
  {
    pos_inhterm = bm@target_var[!(bm@target_var %in% unlist(strsplit(c(arule, irule), '&')))]
  } else
  {
    pos_inhterm = bm@target_var[!(bm@target_var %in% unlist(strsplit(c(arule, irule), '&')) | bm@target_var %in% bm@target_var[index])]
  }

  #pos_inhterm = pos_inhterm[!is.na(pos_inhterm)]
  addlist_irule = list()
  if(length(pos_inhterm)!=0)
  {
    if(irule[1] == '0')
    {
      for(i in 1:length(pos_inhterm))
      {
        addlist_irule = c(addlist_irule, list(pos_inhterm[i]))
      }
    } else
    {
      for(i in 1:length(pos_inhterm))
      {
        addlist_irule = c(addlist_irule, list(c(irule, pos_inhterm[i])))
      }
    }
    
    if(and_bool)
    {
      #Addition of inh rule. (double)
      if(irule[1] != '0')
      {
        for(i in 1:length(pos_inhterm))
        {
          for(j in 1:length(irule))
          {
            if(grepl('&', irule[j]))
            {
              next
            } else
            {
              tmp = sprintf('%s&%s', pos_inhterm[i], irule[j])
              addlist_irule = c(addlist_irule, list(c(irule[-j], tmp)))
            }
          }
        }
      }
    }
  }
  
  #Generate a list of modified models.
  bmodel_dellist = list()
  bmodel_addlist = list()
  
  if(length(dellist_arule) != 0)
  {
    for(i in 1:length(dellist_arule))
    {
      tmp_bm = bm
      tmp_bm@rule_act[[index]] = dellist_arule[[i]]
      if(encoded)
      {
        tmp_bm = list(compress_bmodel(tmp_bm, model_encoding))
        stopifnot(length(model_encoding[[2]])==max(round(tmp_bm[[1]])))
      }
      bmodel_dellist = c(bmodel_dellist, tmp_bm)
    }
  }

  if(length(dellist_irule) != 0)
  {
    for(i in 1:length(dellist_irule))
    {
      tmp_bm = bm
      tmp_bm@rule_inh[[index]] = dellist_irule[[i]]
      if(encoded)
      {
        tmp_bm = list(compress_bmodel(tmp_bm, model_encoding))
        stopifnot(length(model_encoding[[2]])==max(round(tmp_bm[[1]])))
      }
      bmodel_dellist = c(bmodel_dellist, tmp_bm)
    }
  }
  
  if(length(addlist_arule) != 0)
  {
    for(i in 1:length(addlist_arule))
    {
      tmp_bm = bm
      tmp_bm@rule_act[[index]] = addlist_arule[[i]]
      if(encoded)
      {
        tmp_bm = list(compress_bmodel(tmp_bm, model_encoding))
        stopifnot(length(model_encoding[[2]])==max(round(tmp_bm[[1]])))
      }
      bmodel_addlist = c(bmodel_addlist, tmp_bm)
    }
  }
  
  if(length(addlist_irule) != 0)
  {
    for(i in 1:length(addlist_irule))
    {
      tmp_bm = bm
      tmp_bm@rule_inh[[index]] = addlist_irule[[i]]
      if(encoded)
      {
        tmp_bm = list(compress_bmodel(tmp_bm, model_encoding))
        stopifnot(length(model_encoding[[2]])==max(round(tmp_bm[[1]])))
      }
      bmodel_addlist = c(bmodel_addlist, tmp_bm)
    }
  }
  
  return(list(dellist=bmodel_dellist,
              addlist=bmodel_addlist))
}

#' @title Add extra genes to a Boolean model
#' 
#' @description
#' This function adds extra genes to a Boolean model. Return a list of BoolModel object and an initial state.
#' 
#' @param in_gene character vector. Genes to be added into the model.
#' @param in_model data frame or BoolModel object. If it is a data frame, it must have 2 columns, which are targets and update functions.
#' 
#' @export
grow_bmodel = function(in_gene, in_model)
{
  if(class(in_model)=='BoolModel')
  {
    in_model = bm_to_df(in_model)
  }
  
  #Generate a new data frame to be added into the model df.
  empty_func = '(0) &! (0)'
  in_row = data.frame(in_gene, empty_func)
  colnames(in_row) = c('targets', 'factors')
  
  #Modify the model df.
  out_model = rbind(in_model, in_row)
  
  #Convert the model df into BoolModel object.
  out_model = df_to_bm(out_model)

  return(out_model)
}