#' @title Generate random act and inh rule for a single gene
#' 
#' @description
#' This function generates one random Boolean rule (both act and inh) per run. Return a list of two vectors.
#' Note that this method will not give empty rule, i.e. 0 term in both act and inh rules.
#' 
#' @param x character vector. A vector of all single terms to be used.
#' @param np integer. Number of gene variables in a rule. NOT max_varperrule here.
#' @param tar_ind numerical. Indicate which gene is the rule for. Used in preventing self-loop.
#' @param and_bool logical. Indicates whether to include AND terms or not. 
#' @param self_loop logical. Indicates whether to allow self_loop. Default to F.
#' @param rule_type character. Types of rules. Defaults to random.
gen_singlerule = function(x, np, tar_ind, and_bool, self_loop=F, rule_type='random') 
{
  #Convert x to the form of 'v1s'.
  if(!(all(grepl('v[0-9]+s', x))))
  {
    if(any(grepl('v[0-9]+s', x)))
    {
      stop('The name of some genes resemble internally used variables. Respecify names.')
    } else
    {
      x = paste('v', seq(1, length(x)), 's', sep='')
    }
  }
  
  if(rule_type=='full')
  {
    if(self_loop)
    {
      rnum_act = np
      rnum_inh = 0
    } else
    {
      rnum_act = np - 1
      rnum_inh = 0
    }
  } else if(rule_type=='minimal')
  {
    rnum_act = 1
    rnum_inh = 0
  } else if(rule_type=='random')
  {
    rnum_act = floor(runif(1, 0, np+1)) #get a random number from 0 to np for act_rule.
    
    if(rnum_act==0)
    {
      #rnum_inh = floor(runif(1, 1, np-rnum_act+1)) #must have at least one term if rnum_act is empty.
      rnum_inh = np
    } else {
      #rnum_inh = floor(runif(1, 0, np-rnum_act+1)) #get a random number from 0 to np-rnum_act for inh_rule.
      rnum_inh = np - rnum_act
    }
    
    if(rnum_act + rnum_inh > length(x)) #cannot use more terms than available.
    {
      stop('Error in generating rules for random models.')
    }
  }
  
  #Get the activation rule.
  if(rnum_act != 0)
  {
    if(!self_loop)
    {
      ract_prob = rep(1/(length(x)-1), length(x))
      ract_prob[tar_ind] = 0 #set probability to 0, to prevent self-loop.
      tmp_ract = sample(x, rnum_act, prob=ract_prob) #get single variables first.
    } else
    {
      tmp_ract = sample(x, rnum_act) #get single variables first.
    }
    
    rule_act = c()
    
    if(and_bool)
    {
      #Determine which single variables to combine into double variables.
      ra_ind = as.logical(replicate(rnum_act,round(runif(1))))
      while(sum(ra_ind)%%2==1)
      {
        ra_ind = as.logical(replicate(rnum_act,round(runif(1))))
      }
      
      #Add double variables.
      if(length(tmp_ract[ra_ind]) != 0)
      {
        for(i in 1:length(tmp_ract[ra_ind]))
        {
          if(i %% 2 == 1) #only do it when i is odd.
          {
            rule_act = c(rule_act, paste(c(tmp_ract[ra_ind][i], tmp_ract[ra_ind][i+1]), collapse='&'))
          }
        }
      }
      
      #Add single variables.
      rule_act = c(tmp_ract[!ra_ind], rule_act)
    } else
    {
      #Add single variables.
      rule_act = tmp_ract
    }
  } else
  {
    rule_act = '0'
  }
  
  if(rnum_inh != 0)
  {
    #Get the inhibition rule.
    if(!self_loop)
    {
      if(rule_type=='full')
      {
        rinh_prob = rep(1/(length(x)-1), length(x))
        rinh_prob[tar_ind] = 0 #set probability to 0, to prevent self-loop.
        tmp_rinh = sample(x, rnum_inh, prob=rinh_prob) #get single variables first.
      } else
      {
        intr_x = x[!(x %in% unlist(strsplit(rule_act, '&')))]
        rinh_prob = rep(1/(length(intr_x)-1), length(intr_x))
        rinh_prob[which(intr_x %in% x[tar_ind])] = 0 #set probability to 0, to prevent self-loop.
        tmp_rinh = sample(intr_x, rnum_inh, prob=rinh_prob) #get single variables first.
      }
    } else
    {
      tmp_rinh = sample(x[!(x %in% unlist(strsplit(rule_act, '&')))], rnum_inh) #exclude terms already used in rule_act.
    }

    if(and_bool)
    {
      #Determine which single variables to combine into double variables.
      ri_ind = as.logical(replicate(rnum_inh, round(runif(1))))
      while(sum(ri_ind)%%2==1)
      {
        ri_ind = as.logical(replicate(rnum_inh, round(runif(1))))
      }
      
      #Add double variables.
      rule_inh = c()
      if(length(tmp_rinh[ri_ind]) != 0)
      {
        for(i in 1:length(tmp_rinh[ri_ind]))
        {
          if(i %% 2 == 1) #only do it when i is odd.
          {
            rule_inh = c(rule_inh, paste(c(tmp_rinh[ri_ind][i], tmp_rinh[ri_ind][i+1]), collapse='&'))
          }
        }
      }
      
      #Add single variables.
      rule_inh = c(tmp_rinh[!ri_ind], rule_inh)
    } else
    {
      rule_inh = tmp_rinh
    }
  } else
  {
    rule_inh = '0'
  }
  
  return(list(rule_act, rule_inh))
}

#' @title Generate a random Boolean model
#' 
#' @description
#' This function generates a random Boolean model. Returns an S4 BoolModel object.
#' Note that this method will not give empty rule, i.e. 0 term in both act and inh rules.
#' 
#' @param var character vector. A vector of single genes/variables to be used in the model.
#' @param exponent integer. The exponent of power law distribution. Default to 3.
#' @param and_bool logical. Indicates whether to include AND terms or not.
#' @param self_loop logical. Indicates whether to allow self_loop. Default to F.
#' @param mvar integer. Maximum number of variables in act or inh rule. Default to length(var).
#' @param model_type character. Specifies the type of model generated.
#' 
#' @details
#' The number of terms in a function for a gene is modelled by power-law distribution.
gen_one_rmodel = function(var, exponent=3, and_bool, self_loop=F, mvar=length(var), model_type='random')
{
  if(mvar > length(var))
  {
    mvar=length(var)
  }
  
  #(1) By using power law distribution, estimate the in-degree genetic partner for each gene. 
  #The minimum in degree is set to 2. 'v1s' and 'v1s&v2s' are currently considered as 2 different partners.
  if(model_type=='full')
  {
    mvar = length(var)
    num_partner = rep(mvar, length(var))
  } else if(model_type=='minimal')
  {
    mvar = 1
    num_partner = rep(mvar, length(var))
  } else if(model_type=='random')
  {
    num_partner = poweRlaw::rpldis(length(var), xmin=2, alpha=exponent) #xmin = the minimum value of resulting random integer, alpha = scaling factor of the distribution. According to literature, for gene network, this should be 3. (or 2<alpha<3)
    num_partner[num_partner > mvar] = mvar #power law distribution can gives very high number, therefore must cap it.
  }
  
  arule = sapply(1:length(num_partner), function(x) gen_singlerule(var, num_partner[x], x, and_bool, self_loop, rule_type = model_type))
  arule = apply(arule, 1, c) #arule is a list of lists, with list[[1]] = all act rules, list[[2]] = all inh rules.
  
  bmodel = BoolModel(target=var, target_var=paste('v',seq(1,length(var)),'s', sep=''), rule_act=arule[[1]], rule_inh=arule[[2]])
  
  return(bmodel)
}

#' @title Generate two random DAG Boolean models with a specified number of steps apart
#' 
#' @description
#' This function generates a random DAG Boolean model, then get another random DAG Boolean model that is a specified number of steps apart by adding and/or removing genes. 
#' Difficult to generate completely directed graph with a specified number of steps apart.
#' 
#' @param var character vector. A vector of single genes/variables to be used in the model.
#' @param steps integer. Number of steps apart between the two models. If steps=0, give completely random starting model.
#' @param mvar integer. Maximum number of variables in act or inh rule. Default to length(var).
#' @param in_amat matrix. Provide adjacency matrix of first model.
#' @param acyclic logical. Whether to restrict the model to being acyclic or not. Defaults to TRUE.
#' 
#' @export
gen_two_rmodel_dag = function(var, steps, mvar=length(var), in_amat=NULL, acyclic=T)
{ 
  if(!requireNamespace('bnlearn', quietly = TRUE)) {
    stop('Requires bnlearn package to run this function.')
  } else
  {
    if(steps==0)
    {
      start_model = bnlearn::random.graph(var, method = "melancon", max.degree=mvar)
      end_model = bnlearn::random.graph(var, method = "melancon", max.degree=mvar)
    } else
    {
      if(is.null(in_amat))
      {
        #(3) Get first model.
        start_model = bnlearn::random.graph(var, method = "melancon", max.degree=mvar)
      } else
      {
        start_model = bnlearn::empty.graph(var)
        if(acyclic)
        {
          bnlearn::amat(start_model) = in_amat
        } else
        {
          bnlearn::amat(start_model, ignore.cycles=T) = in_amat
        }
      }
      
      cur_ite = 1
      max_ite = 100
      out_break = F
      while(cur_ite < max_ite & !out_break)
      {
        #(4) Get second model of a specified number of steps away.
        end_model = start_model
        memory_ind = c()
        #memory_ind = apply(which(t(amat(end_model))==1, arr.ind = T), 1, function(x) paste(x, collapse=',')) #no undirected arcs.
        #memory_ind = c(memory_ind, paste(1:nrow(tmp), 1:ncol(tmp), sep=',')) #no self_loops.
        #start_time = proc.time()
        while(TRUE)
        {
          tryCatch({
            #used_time = proc.time() - start_time
            #used_time = as.numeric(used_time)[3]
            
            tmp_graph = bnlearn::empty.graph(var)
            tmp = bnlearn::amat(end_model)
            
            x_ind = sample(1:nrow(tmp),1)
            y_ind = sample(1:ncol(tmp),1)
            
            while(paste(x_ind, y_ind, collapse=',') %in% memory_ind)
            {
              x_ind = sample(1:nrow(tmp),1)
              y_ind = sample(1:ncol(tmp),1)
            }
            
            memory_ind = c(memory_ind, paste(x_ind, y_ind, collapse=','))
            if(length(memory_ind) == nrow(tmp)*ncol(tmp)) #break loop, get a new starting model and try again.
            {
              #cat(sprintf('For step %s, no other possible model exists.', steps))
              break
            }
            
            #Change one step at a time.
            if(tmp[x_ind,y_ind]==0)
            {
              tmp[x_ind,y_ind] = 1
            } else
            {
              tmp[x_ind,y_ind] = 0
            }
            
            if(acyclic)
            {
              bnlearn::amat(tmp_graph) = tmp
            } else
            {
              bnlearn::amat(tmp_graph, ignore.cycles=T) = tmp
            }
            
            end_model = tmp_graph
            
            if(sum(abs(bnlearn::amat(start_model)-bnlearn::amat(tmp_graph))) == steps)
            {
              if(nrow(bnlearn::undirected.arcs(tmp_graph))==0) #cannot have undirected arcs in the models.
              {
                out_break = T
                break
              }
            }
          }, error=function(e){},
          finally={})
        }
        
        cur_ite = cur_ite + 1
      }
      
      if(cur_ite == max_ite)
      {
        stop(sprintf('For step %s, no other possible model exists.', steps))
      }
    }
    
    return(list(start_model, end_model))
  }
}

#' @title Generate two random Boolean models with a specified number of steps apart
#' 
#' @description
#' This function generates a random Boolean model, then get another random Boolean model that is a specified number of steps apart by adding and/or removing genes. 
#' Returns a list of two S4 BoolModel objects.
#' 
#' @param var character vector. A vector of single genes/variables to be used in the model.
#' @param steps integer. Number of steps apart between the two models. If steps=0, give completely random starting model.
#' @param mvar integer. Maximum number of variables in act or inh rule. Default to length(var).
#' @param and_bool logical. Indicates whether to include AND terms or not.
#' @param in_bmodel BoolModel object. The starting model supplied.
#' @param self_loop logical. Indicates whether to allow self_loop. Default to F.
#' 
#' @details
#' The number of terms in a function for a gene is modelled by power-law distribution.
#' 
#' @export
gen_two_rmodel = function(var, steps, mvar=length(var), and_bool, in_bmodel=NULL, self_loop=F)
{
  if(mvar > length(var))
  {
    mvar=length(var)
  }
  
  if(is.null(in_bmodel))
  {
    bmodel1 = gen_one_rmodel(var, mvar, and_bool, self_loop)
  } else
  {
    if(class(in_bmodel)!='BoolModel')
    {
      stop('in_bmodel: Please give BoolModel object.')
    }
    
    bmodel1 = in_bmodel
    stopifnot(length(in_bmodel@target) == length(var))
  }
  
  if(steps==0)
  {
    bmodel2 = gen_one_rmodel(var, mvar, and_bool, self_loop)
  } else
  {
    cur_model = bmodel1
    bmodel2 = bmodel1
    ite = 1
    while(length(unlist(model_setdiff(cur_model, bmodel1, mvar)))<steps)
    {
      #       if(length(unlist(model_setdiff(cur_model, bmodel1, mvar))) > steps)
      #       {
      #         browser()
      #         stop('Error in code.')
      #       }
      
      all_rule = lapply(1:length(cur_model@target), function(x) c(cur_model@rule_act[[x]], cur_model@rule_inh[[x]]))
      all_rule = lapply(all_rule, function(x) x[x!='0'])
      all_rule = lapply(all_rule, function(x) unlist(strsplit(x, '&')))
      
      #Summarise info in bmodel.
      all_len = sapply(all_rule, length) #get number of terms (act and inh) for each target gene
      
      avail_add_ind = which(all_len < mvar) #any rule with less genes than mvar can take more genes.
      avail_del_ind = which(all_len > 1) #any rule with more genes than 1 can have fewer genes.
      
      if(length(avail_add_ind) == 0 & length(avail_del_ind) == 0)
      {
        stop('Model is fully specified. No term can be modified.')
        break
      }
      
      #picking choice for modification.
      if(length(avail_add_ind) == 0)
      {
        mod_choice = 'del'
      } else if(length(avail_del_ind) == 0)
      {
        mod_choice = 'add'
      } else
      {
        mod_choice = sample(c('add', 'del'), 1)
      }
      
      #Modifying models.
      if(mod_choice == 'add')
      {
        if(length(avail_add_ind) == 1)
        {
          rule_aind = avail_add_ind
        } else
        {
          rule_aind = sample(avail_add_ind, 1)
        }

        all_models = minmod_model(cur_model, rule_aind, sep_list=T, and_bool=and_bool, self_loop=self_loop)$addlist
        samp_ind = sample(1:length(all_models), 1)[[1]]
        next_model = all_models[[samp_ind]]
        all_models = all_models[-samp_ind]
        while(length(unlist(model_setdiff(cur_model, next_model, mvar)))!=1)
        {
          samp_ind = sample(1:length(all_models), 1)[[1]]
          next_model = all_models[[samp_ind]]
          all_models = all_models[-samp_ind]
        }
      } else if(mod_choice == 'del')
      {
        if(length(avail_del_ind) == 1)
        {
          rule_dind = avail_del_ind
        } else
        {
          rule_dind = sample(avail_del_ind, 1)
        }
        
        all_models = minmod_model(cur_model, rule_aind, sep_list=T, and_bool=and_bool, self_loop=self_loop)$dellist
        samp_ind = sample(1:length(all_models), 1)[[1]]
        next_model = all_models[[samp_ind]]
        all_models = all_models[-samp_ind]
        while(length(unlist(model_setdiff(cur_model, next_model, mvar)))!=1)
        {
          samp_ind = sample(1:length(all_models), 1)[[1]]
          next_model = all_models[[samp_ind]]
          all_models = all_models[-samp_ind]
        }
      }
      stopifnot(length(unlist(model_setdiff(cur_model, next_model, mvar)))==1)
      cur_model = next_model
    }
    
    bmodel2 = cur_model
    stopifnot(length(unlist(model_setdiff(bmodel2, bmodel1, mvar)))==steps)
  }
  
  return(list(bmodel1, bmodel2))
}