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

#' @title Raw single cell qRT-PCR expression data from Wilson et al.
#' 
#' @description
#' A raw single cell expression data obtained from multiple cell types.
#' 
#' @format
#' A data frame with 1626 rows and 44 columns. 
#' 
#' Rows: each row consists of raw expression values from 1 cell.
#' Columns: each column is for 1 gene/variable.
#' 
#' @docType data
#' @name wilson_raw_data
#' @usage data(wilson_raw_data)
NULL

#' @title Raw single cell RNAseq expression data from Wilson et al.
#' 
#' @description
#' A raw single cell expression data obtained from multiple cell types.
#' 
#' @format
#' A data frame with 96 rows and 38498 columns. 
#' 
#' Rows: each row consists of raw expression values from 1 cell.
#' Columns: each column is for 1 gene/variable.
#' 
#' @docType data
#' @name wilson_raw_rnaseq
#' @usage data(wilson_raw_rnaseq)
NULL

#' @title Example Boolean Model used in the vignette
#' 
#' @description
#' A Boolean model used in the examples of the vignette.
#' 
#' @format
#' Each Boolean model is a BoolModel object.
#' 
#' @docType data
#' @name emodel1
#' @usage data(example_models)
NULL

#' @title Example Boolean Model used in the vignette
#' 
#' @description
#' A Boolean model used in the examples of the vignette.
#' 
#' @format
#' Each Boolean model is a BoolModel object.
#' 
#' @docType data
#' @name emodel2
#' @usage data(example_models)
NULL

#' @title Example Boolean Model used in the vignette
#' 
#' @description
#' A Boolean model used in the examples of the vignette.
#' 
#' @format
#' Each Boolean model is a BoolModel object.
#' 
#' @docType data
#' @name emodel3
#' @usage data(example_models)
NULL
