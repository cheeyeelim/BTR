#' @title Output a Boolean Model into Cytoscape readable format
#' 
#' @description
#' This function outputs a Boolean Model in a format that is readable by Cytoscape. Return invisibly the edges (with edge attributes) and node attributes. (i.e. list of 2 dfs)
#' 
#' @param bmodel S4 BoolModel object.
#' @param filepath character vector. Specify path (AND NOT file name). Default to current working directory, i.e. getwd(). Set to NULL to disable file output.
#' 
#' @export
outcyto_model = function(bmodel, filepath=getwd())
{
  edge_vec = character() #setup output vector.
  
  arule = bmodel@rule_act
  irule = bmodel@rule_inh
  
  #Replace the gene names back into the rules.
  for(i in 1:length(bmodel@target_var))
  {
    for(j in 1:length(bmodel@target_var))
    {
      arule[[j]] = gsub(bmodel@target_var[i], bmodel@target[i], arule[[j]])
      irule[[j]] = gsub(bmodel@target_var[i], bmodel@target[i], irule[[j]])
    }
  }
  
  #Split rules into single or double terms.
  s_arule = sapply(arule, function(x) x[!grepl('&', x)])
  d_arule = sapply(arule, function(x) x[grepl('&', x)])
  
  s_irule = sapply(irule, function(x) x[!grepl('&', x)])
  d_irule = sapply(irule, function(x) x[grepl('&', x)])
  
  #Put '0' for empty rules.
  s_arule[sapply(s_arule, function(x) length(x)==0)] = '0'
  d_arule[sapply(d_arule, function(x) length(x)==0)] = '0'
  
  s_irule[sapply(s_irule, function(x) length(x)==0)] = '0'
  d_irule[sapply(d_irule, function(x) length(x)==0)] = '0'
  
  #Start writing rules out for single terms.
  for(i in 1:length(bmodel@target))
  {
    edge_vec = c(edge_vec, paste(s_arule[[i]], 'activates', bmodel@target[i], sep=','))
    edge_vec = c(edge_vec, paste(s_irule[[i]], 'inhibits', bmodel@target[i], sep=','))
  }
  
  #Write rules out for double terms.
  tmp_and = paste('and', seq(1, sum(grepl('&', unlist(d_arule)))+sum(grepl('&', unlist(d_irule)))), sep='_')
  ind = 1
  for(i in 1:length(bmodel@target))
  {
    #Write arule.
    if(d_arule[[i]][1]!='0')
    {
      tmp_rule = strsplit(d_arule[[i]], '&')
      
      for(j in 1:length(tmp_rule))
      {
        #Join coproteins with ANDs.
        edge_vec = c(edge_vec, paste(tmp_rule[[j]], 'activates', tmp_and[ind], sep=','))
        
        #Join ANDs with target genes.
        edge_vec = c(edge_vec, paste(tmp_and[ind], 'activates', bmodel@target[i], sep=','))
        
        ind = ind + 1
      }
    }
    
    #Write irule.
    if(d_irule[[i]][1]!='0')
    {
      tmp_rule = strsplit(d_irule[[i]], '&')
      
      for(j in 1:length(tmp_rule))
      {
        #Join coproteins with ANDs.
        edge_vec = c(edge_vec, paste(tmp_rule[[j]], 'activates', tmp_and[ind], sep=',')) #note that this should be activates.
        
        #Join ANDs with target genes.
        edge_vec = c(edge_vec, paste(tmp_and[ind], 'inhibits', bmodel@target[i], sep=','))
        
        ind = ind + 1
      } 
    }
  }
  
  #Remove anything with 0 in the first and final terms.
  edge_vec = edge_vec[!(grepl('^0,', edge_vec))]
  edge_vec = edge_vec[!(grepl(',0$', edge_vec))]
  
  #Convert the vector into a matrix.
  edge_df = data.frame(do.call(rbind, strsplit(edge_vec, ',')))
  colnames(edge_df) = c('start_node', 'interaction', 'end_node')
  
  #Generating node attributes. (to distinguish gene nodes from AND nodes)
  node_vec = c(bmodel@target, tmp_and)
  node_df = data.frame(node_names=node_vec, node_types=ifelse(grepl('and', node_vec), 'ands', 'genes'))
  
  #Output into files.
  if(!is.null(filepath))
  {
    write.csv(edge_df, file=paste(filepath, '/cytoscape_edges.csv', sep=''), quote=F)
    write.csv(node_df, file=paste(filepath, '/cytoscape_nodes.csv', sep=''), quote=F)
  }
  
  invisible(list(edge_df, node_df))
}

#' @title Output a Boolean Model into Genysis readable format
#' 
#' @description
#' This function outputs a Boolean Model in a format that is readable by Genysis. Return invisibly the formatted vector.
#' 
#' @param bmodel S4 BoolModel object.
#' @param filepath character vector. Specify path (AND NOT file name). Default to current working directory, i.e. getwd(). Set to NULL to disable file output.
#' 
#' @export
outgenysis_model = function(bmodel, filepath=getwd())
{
  gene = bmodel@target
  arule = bmodel@rule_act
  irule = bmodel@rule_inh
  
  #Replace the gene names back into update functions.
  for(i in 1:length(bmodel@target_var))
  {
    arule = lapply(arule, function(x) gsub(bmodel@target_var[i], bmodel@target[i], x))
    irule = lapply(irule, function(x) gsub(bmodel@target_var[i], bmodel@target[i], x))
  }
  
  stopifnot(all(!grepl('v[0-9]+s', unlist(strsplit(unlist(arule), '&'))) | 
                  grepl('^0$', unlist(strsplit(unlist(arule), '&'))))) #check if all gene name conversions is successful.
  stopifnot(all(!grepl('v[0-9]+s', unlist(strsplit(unlist(irule), '&'))) | 
                  grepl('^0$', unlist(strsplit(unlist(irule), '&'))))) #check if all gene name conversions is successful.
  
  #Match gene with rule.
  out_vec = c()
  for(i in 1:length(gene))
  {
    if(arule[[i]][1]!='0')
    {
      out_vec = c(out_vec, paste(arule[[i]], gene[i], sep=' -> '))
    }
    
    if(irule[[i]][1]!='0')
    {
      out_vec = c(out_vec, paste(irule[[i]], gene[i], sep=' -| '))
    }
  }
  
  #Output into files.
  if(!is.null(filepath))
  {
    exist_files = list.files(path=filepath, pattern='^genysis_input_[0-9]+.txt')
    filename = sprintf('genysis_input_%s.txt', length(exist_files)+1)
    
    fcon = file(filename, open='w')
    writeLines(out_vec, fcon)
    close(fcon)
  }
  
  invisible(out_vec)
}

#' @title Generate state transition graph readable by Cytoscapes
#' 
#' @description
#' This function generates a state transition graph using a Boolean model and its state space. Each node represent a state. All nodes in this graph is linked by an edge only if the 2 states have different value in only 1 gene. The output is readable by Cytoscape.
#' 
#' @param mstate data frame. It should be a state(row) x gene(column) df.
#' @param bmodel S4 BoolModel object.
#' @param directed logical. Indicates whether to make directed edges or not. Default to FALSE.
#' @param record.both logical. Indicates whether to also record nodes that have no edges. Default to FALSE.
#' @param filepath character vector. Specify path (AND NOT file name). Default to current working directory, i.e. getwd(). Set to NULL to disable file output.
#' 
#' @export
outcyto_stategraph = function(mstate, bmodel, directed=F, record.both=F, filepath=getwd())
{
  cat(sprintf('Generating state transition network...\n'))
  
  adj_cells = list()
  nonadj_cells = list()
  for(i in 1:nrow(mstate)) #Pick each row.
  {
    cat(sprintf('\rComputing state difference using row %s...', rownames(mstate)[i]))
    
    ind = integer()
    non_ind = integer()
    for(j in 1:nrow(mstate))
    {
      gene_diff = which(mstate[i,]!=mstate[j,])
      
      if(length(gene_diff) == 1) #Only those with one state difference.
      {
        if(directed) #Generate directed network.
        {
          ans = eval_bool(bmodel, unlist(mstate[i,]), gene_diff)
          
          if(ans == mstate[j,gene_diff]) 
          {
            ind = c(ind,j) #Record the state step that is possible.
          } else
          {
            non_ind = c(non_ind, j) #Not possible state step.
          }
        } else #Generate indirected network.
        {
          ind = c(ind,j)
        }
      }
    }
    adj_cells[i] = list(ind)
    nonadj_cells[i] = list(non_ind)
  }
  cat('.\n')
  
  net_all = data.frame(start_node=character(1), interaction=character(1), 
                       end_node=character(1), stringsAsFactors=F) #Both 1 and stringsAsFactors are essential, but 1 will give an empty first row.
  for(i in 1:nrow(mstate))
  {
    if(length(adj_cells[[i]] != 0)) #To exclude cell that does not have partner.
    {
      for(j in 1:length(adj_cells[[i]]))
      {
        net_row = as.character(c(rownames(mstate)[i], '(pp)', 
                                 rownames(mstate)[adj_cells[[i]][j]]))
        
        net_all = rbind(net_all, net_row)
      }
    }
  }
  net_all = net_all[2:nrow(net_all),]
  rownames(net_all) = NULL #Reset row names.
  
  cat(sprintf('State graph generated, with %s nodes and %s edges.\n', length(unique(c(net_all[,1], net_all[,3]))), nrow(net_all)))
  
  if(record.both)
  {
    nonnet_all = data.frame(start_node=character(1), interaction=character(1), 
                            end_node=character(1), stringsAsFactors=F) #Both 1 and stringsAsFactors are essential, but 1 will give an empty first row.
    for(i in 1:nrow(mstate))
    {
      if(length(nonadj_cells[[i]] != 0)) #To exclude cell that does not have partner.
      {
        for(j in 1:length(nonadj_cells[[i]]))
        {
          nonnet_row = as.character(c(rownames(mstate)[i], '(pp)', 
                                      rownames(mstate)[nonadj_cells[[i]][j]]))
          
          nonnet_all = rbind(nonnet_all, nonnet_row)
        }
      }
    }
    nonnet_all = nonnet_all[2:nrow(nonnet_all),]
    rownames(nonnet_all) = NULL #Reset row names.
    
    cat(sprintf('Transition not in agreement with model: %s nodes and %s edges.\n', length(unique(c(nonnet_all[,1], nonnet_all[,3]))), nrow(nonnet_all)))
    
    #Output into files.
    if(!is.null(filepath))
    {
      write.csv(net_all, paste(filepath, '/cytoscape_statespace_edges.txt', sep=''), quote=F, row.names=F)
      write.csv(nonnet_all, paste(filepath, '/cytoscape_notstatespace_edges.txt', sep=''), quote=F, row.names=F)
    }
    
    invisible(list(net_all, nonnet_all))
  }
  
  #Output into files.
  if(!is.null(filepath))
  {
    write.csv(net_all, paste(filepath, '/cytoscape_statespace_edges.txt', sep=''), quote=F, row.names=F)
  }
  
  invisible(net_all)
}
