#' @title Initialise raw data
#' 
#' @description
#' This function initialise raw gene expression values in a matrix. Return a list of two matrices: (1) continuous values and (2) binary values.
#' Note that kmeans clustering as binarisation only works well if the data has a bimodal distribution.
#' 
#' @param x matrix. Numeric data of gene expression.
#' @param max_expr character. Specify whether max expression value is the lowest (as in qPCR), or the highest (as in RNAseq and microarray). Option: 'low', 'high'. Default to 'high'.
#' @param uni_thre numerical. Speficy threshold for unimodality test. Default to 0.2.
#' @param scale logical. Whether to scale the data to a range of 0-1. Default to T.
#' 
#' @export
initialise_raw_data = function(x, max_expr='high', uni_thre=0.2, scale=T)
{
  #(1) Convert negative to positive values.
  if(min(x)<0)
  {
    x = x + abs(min(x))
  } else
  {
    x = x - min(x)
  }
  
  stopifnot(min(x)==0)

  if(max_expr=='low')
  {
    #(2) Invert qPCR values. Lowest expression should be close to 0, highest expression should be away from 0.
    x = abs(max(x) - x)
  }

  #(3) Scale values to between 0 and 1.
  #Scaling should be done by each gene, instead of globally. Because clustering is done by gene, and the expression values should reflect that.
  if(scale)
  {
    for(i in 1:ncol(x))
    {
      a = x[,i]
      if(max(a)-min(a) != 0) #when min and max are different. usual case.
      {
        x[,i] = (a-min(a))/(max(a)-min(a))
      } else
      {
        x[,i] = a-min(a)
      }
    }
  }
  
  #(4) Binarise the data by k-means clustering.
  y = matrix(NA, ncol=ncol(x), nrow=nrow(x))
  for(i in 1:ncol(x))
  {
    #Perform unimodality test for each gene.
    uni_test = diptest::dip.test(x[,i])$p.value
    
    if(uni_test > uni_thre)
    {
      cent_measure = median(x[,i])
      if(cent_measure < 0.5)
      {
        y[,i] = rep(0, nrow(y)) #set all values to 0, if the cent_measure is less than 0.5
      } else
      {
        y[,i] = rep(1, nrow(y)) #set all values to 1, if the cent_measure is more than 0.5
      }
    } else
    {
      bin_result = kmeans(x[,i,drop=F], 2) #Note that kmeans do not give 1 to cluster with lower centroid.
      x_bin = unname(bin_result$cluster)
      
      #(5) Decide which cluster is 0 and which is 1.
      if(bin_result$centers[,1][1] < bin_result$centers[,1][2])
      {
        y[,i] = x_bin - 1 #set 1 to 0, 2 to 1.
      } else
      {
        y[,i] = abs(x_bin-2) #set 1 to 1, 2 to 0.
      }
    }
  }
  
  stopifnot(all(dim(x)==dim(y)))
  
  colnames(y) = colnames(x)
  rownames(y) = rownames(x)
  
  return(list(x, y))
}

#' @title Initialise data
#' 
#' @description
#' This function initialises data frame of Boolean state space. Returns initialised data frame.
#' 
#' @param state data frame. It should contain either 0/1 or F/T data of gene expression.
#' @param aslogic logical. It specifies whether to convert the input data into Boolean values. Default to FALSE.
#' 
#' @export
initialise_data = function(state, aslogic=F)
{
  #Store column and row names.
  rn_tmp = rownames(state)
  cn_tmp = colnames(state)
  
  #Convert to logical T/F if aslogic=T.
  if(aslogic)
  {
    #Convert numerical (i.e. 0, 1) to logical type.
    if(nrow(state) != 1) #one row won't work with apply().
    {
      state = apply(state, 2, as.logical)
    } else 
    {
      state = data.frame(matrix(apply(state, 2, as.logical), nrow=1))
    }

    stopifnot(all(!is.na(state))) #Check if contains NA values due to coercion.
    colnames(state) = cn_tmp
  }
  
  #Convert all row and column names to lowercase.
  #rownames(state) = tolower(rn_tmp)
  #colnames(state) = tolower(cn_tmp)
  
  #Order colname names alphabetically.
  if(all(grepl('v[0-9]+s', cn_tmp))) #if variable names are in the form of v1s, v2s, etc.
  {
    state = state[, order(as.numeric(gsub('v([0-9]+)s', '\\1', cn_tmp)))]
  } else if(all(grepl('G[0-9]+', cn_tmp))) #if variable names are in the form of v1s, v2s, etc.
  {
    state = state[, order(as.numeric(gsub('G([0-9]+)', '\\1', cn_tmp)))]
  } else
  {
    state = state[, order(cn_tmp)]
  }
  
  return(state)
}

#' @title Initialise model
#' 
#' @description
#' This function initialises a Boolean model. Returns initialised S4 BoolModel object.
#' Note that the model should only has 1 NOT operator. More than 1 is STRICTLY NOT allowed.
#' 
#' @param init_model data frame of Boolean model. It should contain two columns, targets and functions.
#' 
#' @export
initialise_model = function(init_model)
{
  bmodel = df_to_bm(init_model)
  return(bmodel)
}

#' @title Remove raw data duplicated wrt to the model state
#' 
#' @description
#' This function removes the 'duplicates' in an expression wrt to the model state.
#' 
#' @param dx matrix. Initialised and discretised numeric data of gene expression.
#' @param cx matrix. Initialised and continuous numeric data of gene expression.
#' 
#' @export
unique_raw_data = function(dx, cx)
{
  fac_dx = c()
  for(i in 1:nrow(dx))
  {
    fac_dx = c(fac_dx, paste(dx[i,], collapse=''))
  }
  
  fac_dx = list(as.factor(fac_dx)) #convert into factors.
  
  output = aggregate(cx, fac_dx, mean) #calculate means for each factor group.
  rownames(output) = output[,1]
  output = output[,-1]
  
  return(as.matrix(output))
}
