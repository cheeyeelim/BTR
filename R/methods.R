#' @title Print Boolean Model
#' 
#' @description
#' This method converts the S4 BoolModel object back into a human-readable data frame, with two columns: (1) target genes, (2) Boolean rules.
#' 
#' @param bmodel S4 BoolModel object.
#' @param gene.names logical. Specify whether to write rules in terms of genes or internal variables. Default to FALSE.
#' 
#' @export
printBM = function(bmodel, gene.names=F)
{
  #Prettify update functions.
  act_rule = sapply(bmodel@rule_act, function(x) paste(x, collapse='|'))
  inh_rule = sapply(bmodel@rule_inh, function(x) paste(x, collapse='|'))
  
  if(gene.names)
  {
    if(all(bmodel@target != bmodel@target_var))
    {
      #Replace the gene names back into update functions.
      for(i in 1:length(bmodel@target_var))
      {
        act_rule = gsub(bmodel@target_var[i], bmodel@target[i], act_rule)
        inh_rule = gsub(bmodel@target_var[i], bmodel@target[i], inh_rule)
      }
      
      stopifnot(all(!grepl('v[0-9]+s', unlist(strsplit(unlist(act_rule), '&'))) | 
                      grepl('^0$', unlist(strsplit(unlist(act_rule), '&'))))) #check if all gene name conversions is successful.
      stopifnot(all(!grepl('v[0-9]+s', unlist(strsplit(unlist(inh_rule), '&'))) | 
                      grepl('^0$', unlist(strsplit(unlist(inh_rule), '&'))))) #check if all gene name conversions is successful.
    }
  }
  
  #Add in brackets around & terms.
  act_rule = lapply(act_rule, function(x) gsub('([[:alnum:]]+&[[:alnum:]]+)', '(\\1)', x))
  inh_rule = lapply(inh_rule, function(x) gsub('([[:alnum:]]+&[[:alnum:]]+)', '(\\1)', x))
  
  comb_rule = paste('(', act_rule, ') &! (', inh_rule, ')', sep='')
  
  #Combine all information together.
  out_df = data.frame(gene=bmodel@target, var=bmodel@target_var, func=comb_rule, stringsAsFactors=F)
  
  return(out_df)
}

#' @title Write Boolean Model
#' 
#' @description
#' This method writes the S4 BoolModel object into a CSV file. This method is a wrapper for print.BoolModel. The output is a data frame, with two columns: (1) target genes, (2) Boolean rules.
#' 
#' @param bmodel S4 BoolModel object.
#' @param file file name with path, or a file connection object.
#' @param gene.names logical. Specify whether to write rules in terms of genes or internal variables. Default to FALSE.
#' @param rownames logical. It specifies whether to write row names.
#' 
#' @export
writeBM = function(bmodel, file, gene.names=F, rownames=F)
{
  #Method that writes information from BoolModel object.
  out_df = printBM(bmodel, gene.names)
  
  write.csv(out_df, file=file, quote=F, row.names=F)
}

#' @title Plot Boolean Model
#' 
#' @description
#' This method plots the network underlying Boolean models by using igraph for quick visualisation.
#' Require igraph.
#' 
#' @param bmodel S4 BoolModel object.
#' @param makePlot logical. Whether to make plot or just return the object. Default to T.
#' @param ... Additional parameters to plot.igraph.
#' 
#' @export
plotBM = function(bmodel, makePlot=T, ...)
{
  #Convert to amat.
  am = bm_to_amat(bmodel)
  
  #Convert into a graph.
  g = igraph::graph.adjacency(am, mode='directed', weighted=T)
  
  #Setup edge colour for plotting.
  #Activation = black, inhibition = red
  igraph::E(g)$color = sapply(igraph::E(g)$weight, function(x) ifelse(x==1, 'black', 'red'))
  
  #Setup other colours.
  igraph::V(g)$frame.color = "white"
  igraph::V(g)$color = rgb(255, 165, 0, 200, maxColorValue = 255)
  
  #Setup vertex font size.
  igraph::V(g)$label.cex = 1.5
  
  if(makePlot)
  {
    #Make the plot.
    igraph::plot.igraph(g, layout=igraph::layout_in_circle, ...)
  }
  
  invisible(g)
}

#' @title Convert BoolModel into adjacency matrix
#' 
#' @description
#' This function converts a BoolModel object into an adjacency matrix.
#' 
#' @param x S4 BoolModel object.
#' @param directed logical. Whether to return directed or undirected adjacency matrix. Default to TRUE.
#' 
#' @export
bm_to_amat = function(x, directed=T)
{ 
  #(1) Extract act and inh rules. And removes &.
  target_gene = x@target
  act_rule = lapply(x@rule_act, function(x) unlist(strsplit(x, '&')))
  inh_rule = lapply(x@rule_inh, function(x) unlist(strsplit(x, '&')))
  all_rule = lapply(1:length(target_gene), function(x) c(act_rule[[x]], inh_rule[[x]]))
  all_rule = lapply(all_rule, function(x) x[x!='0']) #removes 0s.
  
  #(2) Convert gene variables into gene names.
  #Replace the gene names back into update functions.
  for(i in 1:length(x@target_var))
  {
    all_rule = lapply(all_rule, function(y) gsub(x@target_var[i], x@target[i], y))
  }
  
  stopifnot(all(!grepl('v[0-9]+s', unlist(strsplit(unlist(all_rule), '&'))) | 
                  grepl('^0$', unlist(strsplit(unlist(all_rule), '&'))))) #check if all gene name conversions is successful.
  
  #(3) Setup output matrix.
  out_mat = matrix(F, nrow=length(target_gene), ncol=length(target_gene))
  colnames(out_mat) = target_gene
  rownames(out_mat) = target_gene
  
  #(4) Fill in the output matrix.
  for(i in 1:length(target_gene))
  {
    out_mat[,i] = target_gene %in% all_rule[[i]]
  }

  if(!directed)
  {
    #(5) Make the output matrix symmetric.
    for(i in 1:length(target_gene))
    {
      tmp = out_mat[i,] | out_mat[,i]
      out_mat[i,] = tmp
      out_mat[,i] = tmp
    }
    
    stopifnot(isSymmetric(out_mat))
    out_mat = out_mat + 0 #convert logical to numeric.
  } else
  {
    out_mat = out_mat + 0 #convert logical to numeric.
    
    #(6) Make inhibitory edges as negative.
    for(i in 1:length(inh_rule))
    {
      for(j in 1:length(inh_rule[[i]]))
      {
        if(inh_rule[[i]][j]!='0')
        {
          from_ind = as.numeric(gsub('v([0-9]+)s', '\\1', inh_rule[[i]][j]))
          to_ind = i
          
          out_mat[from_ind, to_ind] = -1
        }
      }
    }
  }
  
  return(out_mat)
}

#' @title Convert adjacency matrix into BoolModel object
#' 
#' @description
#' This function converts adjacency matrix to BoolModel object. Able to take in adjacency matrix with -1, which encodes for inhibitory interaction.
#' 
#' @param amat matrix. directed adjacency matrix.
#' @param random logical. Randomly assign to either act or inh rules, when the adjacency matrix only has values 0 and 1, but not -1.
#' @export
amat_to_bm = function(amat, random=F)
{
  if(random)
  {
    stopifnot(any(amat!=-1))
    
    #Get all edges.
    edge_df = cbind(which(amat==1, arr.ind = T), rep(NA, nrow(which(amat==1, arr.ind = T))))
    
    #Setup row and column names.
    colnames(edge_df) = c('from', 'to', 'type')
    rownames(edge_df) = NULL
    
    #Randomly assign edges to act or inh rules.
    edge_df[,3] = replicate(nrow(edge_df), sample(c('act', 'inh'), 1))
    
  } else
  {
    #Get all edges.
    act_edge = cbind(which(amat==1, arr.ind = T), rep('act', nrow(which(amat==1, arr.ind = T))))
    inh_edge = cbind(which(amat==-1, arr.ind = T), rep('inh', nrow(which(amat==-1, arr.ind = T))))
    edge_df = rbind(act_edge, inh_edge)
    
    #Setup row and column names.
    colnames(edge_df) = c('from', 'to', 'type')
    rownames(edge_df) = NULL
  }
  
  #Replace indices by variable names.
  var_name = paste('v', seq(1, length(colnames(amat))), 's', sep='')
  names(var_name) = seq(1, nrow(amat))
  edge_df[,1] = var_name[as.character(edge_df[,1])]
  
  #Convert edge data frame to Boolmodel rules.
  act_rule = vector('list', length(var_name))
  inh_rule = vector('list', length(var_name))
  for(i in 1:nrow(edge_df))
  {
    if(edge_df[i,3] == 'act')
    {
      tmp_act = unname(c(act_rule[[as.numeric(edge_df[i,2])]], edge_df[i,1]))
      act_rule[[as.numeric(edge_df[i,2])]] = tmp_act
    } else
    {
      tmp_inh = unname(c(inh_rule[[as.numeric(edge_df[i,2])]], edge_df[i,1]))
      inh_rule[[as.numeric(edge_df[i,2])]] = tmp_inh
    }
  }
  
  #Fill in empty rules.
  act_rule[which(sapply(act_rule, length)==0)] = '0'
  inh_rule[which(sapply(inh_rule, length)==0)] = '0'
  
  #Generate BoolModel object.
  out_model = BoolModel(target=colnames(amat), target_var=var_name, rule_act=act_rule, rule_inh=inh_rule)
  
  return(out_model)
}

#' @title Convert BoolModel object into BoolNet readable data frame
#' 
#' @description
#' This method converts BoolModel object into a data frame, which is readable by BoolNet.
#' 
#' @param bmodel BoolModel object.
#' 
#' @export
bm_to_df = function(bmodel)
{
  out_df = printBM(bmodel, gene.names=T)[,-2]
  colnames(out_df) = c('targets', 'factors')
  
  return(out_df)
}

#' @title Convert a data frame into BoolModel object
#' 
#' @description
#' This method converts a data frame into a BoolModel object.
#' Note that the model should only has 1 NOT operator. More than 1 is STRICTLY NOT allowed.
#' 
#' @param in_df data frame with 2 columns, targets and factors
#' 
#' @export
df_to_bm = function(in_df)
{
  #Setup the initial data frame.
  in_df = as.data.frame(apply(in_df, 2, tolower), stringsAsFactors=F)
  in_df = in_df[order(in_df[,1]),]
  in_df[,2] = gsub('\\s', '', in_df[,2])
  
  #Take out target genes.
  target = in_df[,1]
  
  #Create corresponding simplified variable terms.
  target_var = paste('v', seq(1,length(target)), 's', sep='') #the variable name consists of digits bounded by two characters.
  
  #Separate the activators and the inhibitors.
  rule_list = strsplit(in_df[,2], '&!')
  
  rule_act = character(length(rule_list)) #initialise variable with correct lengths.
  rule_inh = character(length(rule_list))
  for(i in 1:length(rule_list))
  {
    if(length(rule_list[[i]]) == 2)
    {
      rule_act[i] = rule_list[[i]][1]
      rule_inh[i] = rule_list[[i]][2]
    } else if(length(rule_list[[i]]) == 1)
    {
      if(grepl('!', rule_list[[i]][1])) #if the first term in each vector in the list has a !, this means that there is no activator in this vector, but only inhibitor.
      {
        rule_act[i] = '0'
        rule_inh[i] = rule_list[[i]][1]
      } else
      {
        rule_act[i] = rule_list[[i]][1]
        rule_inh[i] = '0'
      }
    } else
    {
      stop('Error in model rule specification.')
    }
  }
  
  rule_act = gsub('[(]0[)]', '0', rule_act) #take away bracketed 0s.
  rule_inh = gsub('[(]0[)]', '0', rule_inh)
  
  rule_inh = gsub('^([!])', '', rule_inh) #take away remaining !.
  
  #Replace genes in rules with simplified variables.
  for(i in 1:length(target_var))
  {
    rule_act = gsub(target[i], target_var[i], rule_act)
    rule_inh = gsub(target[i], target_var[i], rule_inh)
  }
  
  rule_act = unname(sapply(rule_act, extract_term))
  rule_inh = unname(sapply(rule_inh, extract_term))
  
  if(class(rule_act) != 'list')
  {
    rule_act = as.list(rule_act)
  }
  
  if(class(rule_inh) != 'list')
  {
    rule_inh = as.list(rule_inh)
  }
  
  bmodel = BoolModel(target=target, target_var=target_var, rule_act=rule_act, rule_inh=rule_inh)
  
  return(bmodel)
}
