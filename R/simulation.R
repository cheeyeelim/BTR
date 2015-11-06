#' @title Evaluating Boolean rules
#' 
#' @description
#' This function evaluates the Boolean rules (both act and inh) of one gene at a time. Return a logical value for that gene.
#' 
#' @param bmodel S4 BoolModel object.
#' @param val named logical vector. It should contain the values for all genes at that time point. Note that each value in the vector must be named by its corresponding gene name.
#' @param ind integer. It indicates the state of which gene should be computed.
eval_bool = function(bmodel, val, ind)
{
  #Get all genes, and assign with initial values.
  gene_value = val[bmodel@target] #using gene names to get values.
  
  #Make a gene to variable conversion index.
  gene_var = bmodel@target_var
  names(gene_var) =  bmodel@target
  
  names(gene_value) = gene_var[names(gene_value)] #make name conversion

  #Get the activation and inhibition rules and replace gene names with values in rules.
  act_rule = bmodel@rule_act[[ind]]
  inh_rule = bmodel@rule_inh[[ind]]
  
  for(i in names(gene_value))
  {
    act_rule = gsub(i, gene_value[i], act_rule)
    inh_rule = gsub(i, gene_value[i], inh_rule)
  }
  
  #Evaluate the activation and inhibition rules.
  act_ans = eval(parse(text=paste(act_rule,collapse='|'))) #join with OR and evaluate. parse converts string into expression.
  inh_ans = eval(parse(text=paste(inh_rule,collapse='|')))
  
  return(act_ans&!inh_ans)
}

#' @title Simulating Boolean model
#' 
#' @description
#' This function simulates the Boolean model using an initial state. Returns the full asynchronous state space, and point steady states.
#' 
#' @param bmodel S4 BoolModel object.
#' @param istate data frame. It must have been initialised by initialise_data(), and has gene names as column names. Must contain only 1 row.
#' @param steady_bool logical. Specifies whether to return point steady states or not. Default to F.
#' 
#' @export
simulate_model = function(bmodel, istate, steady_bool=F)
{
  stopifnot(nrow(istate)==1) #must have only 1 row. Only 1 initial state should be passed into this function at each iteration.
  
  #Filter the start state, if the model is modified. And convert it into a named vector.
  #this step will filter and ensure that both bmodel and istate have the same gene order.
  istate = istate[bmodel@target]
  istate = unname(unlist(istate))
  
  #Convert the bmodel from R object into R list.
  lmodel = decreate_boolmodel(bmodel)
  
  all_state = rcpp_simulate(lmodel, istate) #note that a list is returned.
  
  #Convert the list into a df.
  full_state = t(as.data.frame(all_state[[1]]))

  #Set row and column names.
  rownames(full_state) = paste('state', seq(1,nrow(full_state)), sep='_')
  colnames(full_state) = bmodel@target
  
  if(steady_bool)
  {
    steady_state = t(as.data.frame(all_state[[2]]))
    if(nrow(steady_state)!=0)
    {
      rownames(steady_state) = paste('state', seq(1,nrow(steady_state)), sep='_')
      colnames(steady_state) = bmodel@target
    }
    
    return(list(full_state, steady_state))
  }
  
  return(full_state)
}

#' @title Decreate Boolean model
#' 
#' @description
#' This function converts a S4 BoolModel object into a list of 4 lists. (for use by Rcpp in simulate_model())
#' 
#' @param bmodel S4 BoolModel object.
decreate_boolmodel = function(bmodel) #boolean_model
{
  out_list = list()
  for(i in 1:length(slotNames(bmodel)))
  {
    tmp_element = slot(bmodel, slotNames(bmodel)[i]) #use slot() to obtain the elements from the specific slot.
    out_list = c(out_list, list(tmp_element))
  }
  
  #Give names to each list within the main list. (For easier access by Rcpp)
  names(out_list) = c('gene','var', 'act_rule', 'inh_rule')
  
  return(out_list)
}
