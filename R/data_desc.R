#' @title HSC Boolean Model from Bonzanni et al.
#' 
#' @description
#' A Boolean model describing HSC in mice. It contains 11 genes. Its steady state is a cyclic loop of 32 states.
#' 
#' @format
#' A data frame with 11 rows and 2 columns. 
#' 
#' Rows: each row consists of 1 gene and its associated Boolean rule.
#' Column 1: target gene
#' Column 2: associated Boolean rule
#' 
#' @docType data
#' @name bon_bmodel
#' @usage data(bon_bmodel)
NULL

#' @title Initial state from Bonzanni et al.
#' 
#' @description
#' An intial state specified in Bonzanni et al. It contains a set of Boolean values for 11 genes.
#' 
#' @format
#' A data frame with 1 row and 11 columns. 
#' 
#' Rows: each row consists of 1 set of Boolean state.
#' Columns: each column is for 1 gene/variable.
#' 
#' @docType data
#' @name bon_istate
#' @usage data(bon_istate)
NULL

#' @title Initial state from Moignard et al.
#' 
#' @description
#' An intial state obtained from data in Moignard et al, determined by taking colMeans over unique rows, and rounding the means to 0-1. 
#' Values for genes that are missing in Moignard et al, but are present in Bonzanni et al, are determined by taking values from the original initial state supplied in Bonzanni et al.
#' It contains a set of Boolean values for 20 genes.
#' 
#' @format
#' A data frame with 1 row and 20 columns. 
#' 
#' Rows: each row consists of 1 set of Boolean state.
#' Columns: each column is for 1 gene/variable.
#' 
#' @docType data
#' @name bon_moig_istate
#' @usage data(bon_moig_istate)
NULL

#' @title Myeloid Boolean Model from Krumsiek et al.
#' 
#' @description
#' A Boolean model describing myeloid development in mice. It contains 11 genes. Its steady states are 4 static attractors.
#' 
#' @format
#' A data frame with 11 rows and 2 columns. 
#' 
#' Rows: each row consists of 1 gene and its associated Boolean rule.
#' Column 1: target gene
#' Column 2: associated Boolean rule
#' 
#' @docType data
#' @name krum_bmodel
#' @usage data(krum_bmodel)
NULL

#' @title Initial state from Krumsiek et al.
#' 
#' @description
#' An intial state specified in Krumsiek et al. It contains a set of Boolean values for 11 genes.
#' 
#' @format
#' A data frame with 1 row and 11 columns. 
#' 
#' Rows: each row consists of 1 set of Boolean state.
#' Columns: each column is for 1 gene/variable.
#' 
#' @docType data
#' @name krum_istate
#' @usage data(krum_istate)
NULL

#' @title Raw single cell qRT-PCR expression data from Moignard et al.
#' 
#' @description
#' A raw single cell expression data obtained from multiple cell types.
#' 
#' @format
#' A data frame with 597 rows and 18 columns. 
#' 
#' Rows: each row consists of raw expression values from 1 cell.
#' Columns: each column is for 1 gene/variable.
#' 
#' @docType data
#' @name moig_raw_data
#' @usage data(moig_raw_data)
NULL

#' @title Discretised single cell qRT-PCR expression data from Moignard et al.
#' 
#' @description
#' A discretised single cell expression data obtained from multiple cell types.
#' 
#' @format
#' A data frame with 597 rows and 18 columns. 
#' 
#' Rows: each row consists of discretised expression values from 1 cell.
#' Columns: each column is for 1 gene/variable.
#' 
#' @docType data
#' @name moig_data
#' @usage data(moig_data)
NULL

#' @title Estimated parameters from Wilson et al. data
#' 
#' @description
#' A list of parameters (based on log normal distribution) estimated from Wilson et al. single-cell qPCR expression data.
#' 
#' @format
#' A list with 4 numeric vectors, all_mu1, all_mu2, all_sig1, all_sig2. Note that each element in the vector is estimated from a single gene.
#' 
#' @docType data
#' @name real_param
#' @usage data(real_param)
NULL