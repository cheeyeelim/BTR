% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialisation.R
\name{initialise_model}
\alias{initialise_model}
\title{Initialise model}
\usage{
initialise_model(init_model)
}
\arguments{
\item{init_model}{data frame of Boolean model. It should contain two columns, targets and functions.}
}
\description{
This function initialises a Boolean model. Returns initialised S4 BoolModel object.
Note that the model should only has 1 NOT operator. More than 1 is STRICTLY NOT allowed.
}

