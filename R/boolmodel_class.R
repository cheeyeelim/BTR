#' @title An S4 class to represent a Boolean Model
#' 
#' @description
#' This class represents Boolean Model in a S4 BoolModel object.
#' 
#' @field target character vector. It should contain gene symbols.
#' @field targer_var character vector. It should contain the internal representation of gene variables, in the form of "v[0-9]+s"
#' @field rule_act list of character vectors. Each element in the list should contain the activating gene variables for a particular target gene.
#' @field rule_inh list of character vectors. Each element in the list should contain the inhibiting gene variables for a particular target gene.
#' 
#' @export BoolModel
#' @exportClass BoolModel
BoolModel = setClass('BoolModel',
                       slots=c(target='character', target_var='character', 
                               rule_act='list', rule_inh='list')
                     )

