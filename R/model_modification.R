#' @title Minimal modification of whole Boolean model
#' 
#' @description
#' This function generates all possible boolean models minimally modified. Returns a lists of 2 lists, deleted models and added models
#' 
#' @param bm S4 BoolModel object.
#' @param index integer. Specifying rule of which gene to modify. If NULL, modifies all rules in the model. Defaults to NULL.
#' @param overlap_gene character vector. Specify which genes are present in both model and data inputs. Only needed when index=NULL.
minmod_model = function(bm, index=NULL, overlap_gene=NULL)
{
  if(is.null(index))
  {
    #Modify only the overlapping genes
    gene_ind = which(bm@target %in% overlap_gene)
    
    del_list = c()
    add_list = c()
    for(ind in gene_ind) #Modify act and inh rules for one gene at a time.
    {
      tmp_list = minmod_internal(bm, ind)
      if(!is.null(tmp_list$dellist))
      {
        del_list = c(del_list, tmp_list$dellist)
      }
      if(!is.null(tmp_list$addlist))
      {
        add_list = c(add_list, tmp_list$addlist)
      }
    }
    out_list = list(del_list=del_list,
                    add_list=add_list)
  } else
  {
    out_list = minmod_internal(bm, index)
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
minmod_internal = function(bm, index)
{
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
    } else if(arule[1] != '0' & irule[1] != '0')
    {
      dellist_arule = c(dellist_arule, list('0'))
    }
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

  #Deletion of inh rule.
  dellist_irule = list()
  if(length(irule) == 1)
  {
    if(grepl('&' , irule))
    {
      tmp = unlist(strsplit(irule[1], '&')) #split the AND term into 2.
      dellist_irule = c(dellist_irule, list(c(tmp[1], irule[-1]), c(tmp[2], irule[-1])))
    } else if(arule[1] != '0' & irule[1] != '0')
    {
      dellist_irule = c(dellist_irule, list('0'))
    }
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
  
  #Addition of act rule. (single)
  pos_actterm = bm@target_var[!bm@target_var %in% unlist(strsplit(c(arule, irule), '&'))]
  addlist_arule = list()
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
  
  #Addition of inh rule. (single)
  pos_inhterm = bm@target_var[!bm@target_var %in% unlist(strsplit(c(arule, irule), '&'))]
  addlist_irule = list()
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
  
  #Generate a list of modified models.
  bmodel_dellist = list()
  bmodel_addlist = list()
  
  if(length(dellist_arule) != 0)
  {
    for(i in 1:length(dellist_arule))
    {
      tmp_bm = bm
      tmp_bm@rule_act[[index]] = dellist_arule[[i]]
      bmodel_dellist = c(bmodel_dellist, tmp_bm)
    }
  }

  if(length(dellist_irule) != 0)
  {
    for(i in 1:length(dellist_irule))
    {
      tmp_bm = bm
      tmp_bm@rule_inh[[index]] = dellist_irule[[i]]
      bmodel_dellist = c(bmodel_dellist, tmp_bm)
    }
  }
  
  if(length(addlist_arule) != 0)
  {
    for(i in 1:length(addlist_arule))
    {
      tmp_bm = bm
      tmp_bm@rule_act[[index]] = addlist_arule[[i]]
      bmodel_addlist = c(bmodel_addlist, tmp_bm)
    }
  }
  
  if(length(addlist_irule) != 0)
  {
    for(i in 1:length(addlist_irule))
    {
      tmp_bm = bm
      tmp_bm@rule_inh[[index]] = addlist_irule[[i]]
      bmodel_addlist = c(bmodel_addlist, tmp_bm)
    }
  }
  
  return(list(dellist=bmodel_dellist,
              addlist=bmodel_addlist))
}

#' @title Add extra genes to a Boolean model
#' 
#' @description
#' This function adds extra genes to a Boolean model. Input model must be in data frame format, output model will be BoolModel object.
#' 
#' @param in_model data frame with 2 columns, which are targets and factors
#' @param in_gene character vector. Genes to be added into the model.
#' 
#' @export
grow_bmodel = function(in_model, in_gene)
{
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